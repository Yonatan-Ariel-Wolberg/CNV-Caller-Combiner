#!/usr/bin/env nextflow

// =============================================================================
// Detailed documentation: dragen_germline_cnv_end_to_end.nf
// Purpose: This script supports CNV calling and/or post-processing tasks in
//   this project, with parameters and I/O conventions aligned to repository workflows.
// Scope: End-to-end run path for DRAGEN.
// Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
// Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
// Operational notes: Keep parameter names synchronized with related workflow/process
//   definitions and update this file alongside upstream interface changes.
// =============================================================================
nextflow.enable.dsl=2
include { RUN_MOSDEPTH; BUILD_CONSENSUS } from './modules/common/modules-mosdepth-target.nf'

include { MERGE_CONSENSUS_COHORT; SPLIT_COHORT_VCF; sanitizeCohortVcf; bcftoolsFlattenCohort } from './modules/consensus/consensus-postprocess.nf'
include { DOWNLOAD_ANNOTSV_DB; runAnnotSV; addACMGAnnotations } from './modules/annotation/annotation.nf'
include { filterRareCnvs; filterVcfByRareRegions } from './modules/filters/allele-frequency-filters.nf'
include { EVALUATE_ACMG_CSV; PLOT_EVALUATION_METRICS } from './modules/evaluation/performance-metric-evaluation.nf'
include { ANALYZE_CNV_QUALITY_METRICS } from './modules/evaluation/quality-metric-evaluation.nf'
include { RUN_CATBOOST_PREDICT } from './modules/ml/modules-predict.nf'
include { FILTER_CANOES_QUALITY } from './modules/filters/quality-filters.nf'
include { FILTER_CLAMMS_QUALITY } from './modules/filters/quality-filters.nf'
include { FILTER_XHMM_QUALITY } from './modules/filters/quality-filters.nf'
include { FILTER_CNVKIT_QUALITY } from './modules/filters/quality-filters.nf'
include { FILTER_GCNV_QUALITY } from './modules/filters/quality-filters.nf'
include { FILTER_DRAGEN_QUALITY } from './modules/filters/quality-filters.nf'
include { NORMALISE_CNV_QUALITY_SCORES } from './modules/normalise/modules-normalise.nf'


/*
 * -------------------------------------------------------------------------------------
 *  DRAGEN GERMLINE CNV END-TO-END ANALYSIS PIPELINE
 *  Consensus BED Generation -> DRAGEN Caller -> Cohort Filtering & ACMG -> Evaluation & Plotting
 * -------------------------------------------------------------------------------------
 */

// 1. Pipeline Modes & Options
params.run_modes             = "Control,Mosdepth" // Modes to run. Can be "Control", "Mosdepth", or "Control,Mosdepth"
params.run_calling           = true  // Execute CNV calling on Illumina Connected Analytics (ICAv2)
params.pass_percent          = 90    // Consensus threshold: regions covered in >= X% of samples

// 2. Reference & Input Paths
params.fasta                 = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.fasta_fai             = "${params.fasta}.fai"
params.standard_bed          = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
params.cram_dir              = "/dataG/ddd/data/ega/4_ddd_truthset_aligned_2026_04_30/ega_truthset/crams"

// 3. ICAv2 Credentials & Configuration (required if params.run_calling = true)
params.projectId             = ""
params.pipelineId            = ""
params.userReference         = "ddduk-pipeline-run"
params.storageSize           = "Medium"
params.icaUploadPath         = "/uploads/crams/"
params.localDownloadPath     = "./raw_calls"

// Analysis Data Codes from params.json
params.referenceAnalysisDataCode      = "ref_tar"
params.targetBedAnalysisDataCode      = "bed"
params.cramAnalysisDataCode           = "crams"
params.craiAnalysisDataCode           = "crais"
params.cramReferenceAnalysisDataCode  = "ref_cram"
params.referenceFileId                = ""
params.cramReferenceFileId            = ""
params.samplesPerNode                 = 8

// 4. Pre-called VCF Fallback (used if params.run_calling = false)
params.precalled_control_vcfs  = ""
params.precalled_mosdepth_vcfs = ""
params.precalculated_consensus_bed = "/home/ywolberg/dragen-germline-enrichment-v4-4-6/multiqc-filter-mosdepth/ddduk_karyotypic_consensus.bed"

// 5. Cohort Filtering & ACMG Annotations
params.total_samples         = 90
params.max_af                = 0.05
params.outdir                = "./results"
params.annotsv_db            = ""

// 6. Evaluator Paths
params.map_file              = "sample_id.map"
params.truth_bed             = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed"

// Define dynamic repository and bin directory paths depending on script location
def repoDir = file("${baseDir}/main.nf").exists() ? "${baseDir}" : (file("${baseDir}/../main.nf").exists() ? "${baseDir}/.." : "${baseDir}")
def binDir = file("${repoDir}/bin").exists() ? "${repoDir}/bin" : "${repoDir}"


/*
 * =====================================================================================
 *  PHASE 1: MOSDEPTH TARGET BED GENERATION
 * =====================================================================================
 */






/*
 * =====================================================================================
 *  PHASE 2: ICAv2 DRAGEN SEQUENTIAL CALLING WORKFLOW
 *  We use maxForks 1 to guarantee that only one run passes through ICAv2 calling, 
 *  downloads outputs, and deletes its inputs/outputs before the next run starts.
 * =====================================================================================
 */

process RUN_ICA_CALLING_WORKFLOW {
    tag "${run_name}"
    maxForks 1
    cpus 1
    memory '8 GB'
    publishDir "${params.outdir}/${run_name}/raw_calls", mode: 'copy'

    input:
    tuple val(run_name), path(bed_file)
    path cram_files

    output:
    tuple val(run_name), path('raw_vcf_outputs/**/*.cnv.vcf.gz'), emit: vcfs

    script:
    def projectId = params.projectId
    def uploadPath = params.icaUploadPath
    def pipelineId = params.pipelineId
    def storageSize = params.storageSize
    def userRefPrefix = params.userReference
    def samplesPerNode = params.samplesPerNode
    
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.craiAnalysisDataCode
    def refCode = params.referenceAnalysisDataCode
    def bedCode = params.targetBedAnalysisDataCode
    def cramRefCode = params.cramReferenceAnalysisDataCode
    def refId = params.referenceFileId
    def cramRefId = params.cramReferenceFileId
    """
    #!/bin/bash
    set -euo pipefail
    
    clean_json() { sed -n '/^{/,\$p'; }

    get_or_upload_file() {
        local fpath=\$1
        local folder=\$2
        local pid=\$3
        local fname=\$(basename "\$fpath")
        
        local list_raw=\$(icav2 projectdata list --project-id "\$pid" --parent-folder "\$folder" --file-name "\$fname" --match-mode EXACT -o json 2>/dev/null || true)
        local existing_id=\$(echo "\$list_raw" | clean_json | jq -r --arg n "\$fname" '.items[]? | select(.name == \$n) | .id' | head -n 1 || true)
        local existing_status=\$(echo "\$list_raw" | clean_json | jq -r --arg n "\$fname" '.items[]? | select(.name == \$n) | .status' | head -n 1 || true)
        
        if [ -n "\$existing_id" ] && [ "\$existing_id" != "null" ]; then
            if [ "\$existing_status" == "PARTIAL" ]; then
                icav2 projectdata delete "\$existing_id" --project-id "\$pid" >/dev/null 2>&1 || true
            else
                local chk=\$(icav2 projectdata get "\$existing_id" --project-id "\$pid" -o json 2>/dev/null || true)
                if ! echo "\$chk" | grep -qE "NotFound|not found|404"; then
                    echo "\$existing_id"
                    return 0
                fi
            fi
        fi
        
        local real_fpath=\$(realpath "\$fpath")
        local upload_raw=\$(timeout 3600 icav2 projectdata upload "\$real_fpath" "\$folder" --project-id "\$pid" -o json 2>&1 || true)
        
        local new_id=\$(echo "\$upload_raw" | clean_json | jq -r '.id // empty' 2>/dev/null || true)
        if [ -z "\$new_id" ] || [ "\$new_id" == "null" ]; then 
            new_id=\$(echo "\$upload_raw" | grep -o 'fil.[a-zA-Z0-9]*' | head -1 || true)
        fi
        echo "\$new_id"
    }

    # 1. Upload Target BED file dynamically
    echo "Uploading target BED file for ${run_name}..."
    BED_ID=\$(get_or_upload_file "${bed_file}" "${uploadPath}" "${projectId}" | tr -d '[:space:]')
    echo "Uploaded BED ID: \$BED_ID"

    # 2. Upload CRAM & CRAI files in parallel (concurrency of 8)
    echo "Uploading CRAM and CRAI files..."
    mkdir -p upload_results
    cram_list=( ${cram_files} )

    upload_task() {
        local file_path=\$1
        local fname=\$(basename "\$file_path")
        local file_id
        file_id=\$(get_or_upload_file "\$file_path" "${uploadPath}" "${projectId}" | tr -d '[:space:]')
        echo "\$fname:\$file_id" > "upload_results/\$fname.txt"
    }
    export -f get_or_upload_file clean_json upload_task

    for f in "\${cram_list[@]}"; do
        upload_task "\$f" &
        while [ \$(jobs -r | wc -l) -ge 8 ]; do
            sleep 2
        done
    done
    wait

    # Compile file IDs
    cram_ids=""
    crai_ids=""
    for f in "\${cram_list[@]}"; do
        fname=\$(basename "\$f")
        res=\$(cat "upload_results/\$fname.txt")
        fid=\$(echo "\$res" | cut -d: -f2)
        if [[ "\$fname" == *.cram ]]; then
            cram_ids="\${cram_ids:+\$cram_ids,}\$fid"
        else
            crai_ids="\${crai_ids:+\$crai_ids,}\$fid"
        fi
    done

    # 3. Start batch CNV calling analysis
    echo "Launching DRAGEN batch CNV calling on ICAv2..."
    user_ref="${userRefPrefix}-${run_name}"

    analysis_response_raw=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --project-id ${projectId} \
        --user-reference "\${user_ref}" \
        --storage-size ${storageSize} \
        --input ${refCode}:${refId} \
        --input ${bedCode}:\$BED_ID \
        --input ${cramRefCode}:${cramRefId} \
        --input ${cramCode}:\$cram_ids \
        --input ${craiCode}:\$crai_ids \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:true \
        --parameters output_format:CRAMv3.1 \
        --parameters enable_duplicate_marking:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:GVCF \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:true \
        --parameters cnv_enable_in_run_pon:true \
        --parameters cnv_enable_gcbias_correction:true \
        --parameters enable_sv:true \
        --parameters qc_coverage_report:cov_report \
        --parameters qc_coverage_filters:"'mapq<20,bq<20'" \
        --parameters qc_detect_contamination:true \
        --parameters enable_variant_annotation:true \
        --parameters samples_per_node:${samplesPerNode} 2>&1 || true)
        
    analysis_id=\$(echo "\$analysis_response_raw" | sed -n '/^{/,\$p' | jq -r ".id" 2>/dev/null || echo "")
    if [ -z "\$analysis_id" ] || [ "\$analysis_id" == "null" ]; then
        echo "Error: Failed to launch analysis batch."
        echo "\$analysis_response_raw"
        exit 1
    fi
    echo "Analysis ID: \$analysis_id"

    # 4. Monitor analysis status until finished
    echo "Monitoring analysis..."
    while true; do
        status_json=\$(icav2 projectanalyses get "\$analysis_id" --project-id "${projectId}" -o json)
        status=\$(echo "\$status_json" | jq -r ".status" | tr -d '[:space:]')
        echo "Current Status: \$status"
        if [[ "\$status" == "SUCCEEDED" ]]; then
            break
        elif [[ "\$status" == "FAILED" || "\$status" == "FAILED_FINAL" || "\$status" == "ABORTED" ]]; then
            echo "Error: ICAv2 Analysis failed."
            exit 1
        fi
        sleep 300
    done

    # 5. Download CNV calling outputs locally
    echo "Downloading calling outputs..."
    mkdir -p raw_vcf_outputs
    output_json=\$(icav2 projectanalyses output "\$analysis_id" --project-id "${projectId}" -o json)
    folder_id=\$(echo "\$output_json" | jq -r ".items[0].data[0].dataId")

    if [ -n "\$folder_id" ] && [ "\$folder_id" != "null" ]; then
        icav2 projectdata download "\$folder_id" "raw_vcf_outputs/" --project-id "${projectId}"
    else
        echo "Error: Output folder ID not found."
        exit 1
    fi

    # 6. Cleanup ICAv2 uploaded files and outputs
    echo "Deleting inputs and outputs off ICAv2 to clear space..."
    for fid in \${cram_ids//,/ }; do
        icav2 projectdata delete "\$fid" --project-id "${projectId}" || true
    done
    for fid in \${crai_ids//,/ }; do
        icav2 projectdata delete "\$fid" --project-id "${projectId}" || true
    done
    icav2 projectdata delete "\$folder_id" --project-id "${projectId}" || true
    icav2 projectdata delete "\$BED_ID" --project-id "${projectId}" || true
    """
}


/*
 * =====================================================================================
 *  PHASE 3: COHORT FILTERING, MERGING & ACMG PRIORITIZATION
 * =====================================================================================
 */



process SURVIVORMergeCohort {
    label 'survivor'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(vcf_list), path(filtered_vcfs)

    output:
    tuple val(run_name), path('cohort_merged.vcf'), emit: merged_vcf

    script:
    """
    SURVIVOR merge vcf_list.txt 1000 1 1 0 0 0 cohort_merged.vcf
    """
}
















/*
 * =====================================================================================
 *  PHASE 4: RIGOROUS EVALUATION & VISUALIZATION
 * =====================================================================================
 */






/*
 * =====================================================================================
 *  WORKFLOW ENTRY POINT
 * =====================================================================================
 */




workflow {
    if (workflow.stubRun) {
        log.info "Validation-only execution mode detected (stub-run); exiting before process execution."
        return
    }

    // 1. Reference Validation & Inputs Setup
    def modes = params.run_modes.tokenize(',')
    
    // Validate that inputs for mosdepth are present if Mosdepth run mode is selected
    if (modes.contains('Mosdepth')) {
        if (!params.bam_dir) error "DRAGEN consensus generation (Mosdepth mode) requires --bam_dir containing BAM files."
    }

    // Initialize the AnnotSV database path channel
    if (params.annotsv_db) {
        annotsv_db_ch = Channel.value(file(params.annotsv_db))
    } else {
        def local_db = file("${repoDir}/annotsv_db/AnnotSV_annotations")
        if (local_db.exists()) {
            annotsv_db_ch = Channel.value(local_db)
        } else {
            DOWNLOAD_ANNOTSV_DB()
            annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir.first()
        }
    }

    // 2. Mosdepth consensus BED generation (only executed if Mosdepth mode is active)
    if (modes.contains('Mosdepth')) {
        if (params.precalculated_consensus_bed) {
            mosdepth_bed_ch = Channel.fromPath(params.precalculated_consensus_bed)
        } else {
            cram_ch = Channel.fromPath("${params.bam_dir}/*.bam")
                .map { bam -> 
                    def bai = file("${bam}.bai")
                    if (!bai.exists()) bai = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
                    if (!bai.exists()) error "Missing BAM index file for ${bam}"
                    return [bam.baseName, bam, bai]
                }
            
            mosdepth_out = RUN_MOSDEPTH(
                cram_ch, 
                file(params.fasta), 
                file(params.fasta_fai), 
                file(params.standard_bed)
            )
            
            BUILD_CONSENSUS(mosdepth_out.usable_beds.collect(), file(params.standard_bed))
            mosdepth_bed_ch = BUILD_CONSENSUS.out.karyotypic_bed
        }
    } else {
        mosdepth_bed_ch = Channel.empty()
    }

    // 3. Construct target_beds_ch pairing run_name with target BED
    control_bed_ch = modes.contains('Control') ? Channel.of([ 'Control', file(params.standard_bed) ]) : Channel.empty()
    mosdepth_bed_entry_ch = modes.contains('Mosdepth') ? mosdepth_bed_ch.map { bed -> [ 'Mosdepth', bed ] } : Channel.empty()
    
    target_beds_ch = control_bed_ch.concat(mosdepth_bed_entry_ch)

    // 4. CNV Variant Calling Phase (Dynamic Calling vs Precalled Fallback)
    if (params.run_calling) {
        called_modes_ch = target_beds_ch.filter { run_name, bed ->
            def existing = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            if (existing && existing.size() > 0) {
                log.info "Bypassing ICA DRAGEN calling for mode ${run_name} because raw calls were found in results directory."
                return false
            }
            return true
        }

        bypassed_vcfs_ch = target_beds_ch.filter { run_name, bed ->
            def existing = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            return (existing && existing.size() > 0)
        }.map { run_name, bed ->
            def files = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            return [ run_name, files ]
        }

        // Fetch BAM files to run DRAGEN sequentially
        cram_files_ch = Channel.fromPath(params.bam_dir + '/*.bam')
            .map { f -> [f.baseName, [f, file(f.toString() + '.bai')]] }
        
        // Execute caller sequentially per run mode (maxForks 1 guarantees serial execution and clean deletes)
        RUN_ICA_CALLING_WORKFLOW(called_modes_ch, cram_files_ch.collect())
        vcfs_ch = RUN_ICA_CALLING_WORKFLOW.out.vcfs.concat(bypassed_vcfs_ch)
    } else {
        // Bypass calling and select pre-called VCF sets
        vcfs_ch = Channel.fromList(modes)
            .map { run_name ->
                def glob_path = (run_name == 'Control') ? params.precalled_control_vcfs : params.precalled_mosdepth_vcfs
                def files = file(glob_path)
                return [ run_name, files ]
            }
    }

    // 5. Cohort Analysis & ACMG Priorities (Runs in parallel across run modes)
    FILTER_DRAGEN_QUALITY(vcfs_ch)

    if (!params.skip_annotation) {
        SURVIVORMergeCohort(FILTER_DRAGEN_QUALITY.out.filtered_data)
        bcftoolsFlattenCohort(SURVIVORMergeCohort.out.merged_vcf)
        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
        
        // Join merged VCF with rare regions
        slice_joined = SURVIVORMergeCohort.out.merged_vcf.join(filterRareCnvs.out.rare_data)
        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
        
        runAnnotSV(filterVcfByRareRegions.out.rare_vcf, annotsv_db_ch)
        
        // Join unannotated rare CSV with AnnotSV TSV output
        acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })

        // 6. CNV Quality Metrics Assessment
        // Group all inputs for the metrics analysis
        metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
            .join(runAnnotSV.out.annotated_tsv)
            .join(FILTER_DRAGEN_QUALITY.out.filtered_data.map { [ it[0], it[2] ] })
            
        // Combine with target BED files to get the correct bed file per run_name
        metrics_joined = metrics_inputs.combine(target_beds_ch, by: 0)
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_joined.map { [ it[0], it[1], it[2], it[3] ] }, // [run_name, rare_cnvs_csv, annotated_tsv, filtered_vcfs]
            metrics_joined.map { it[4] } // bed_file
        )

        // 7. Score-separated Evaluation
        // Expand the list of ACMG CSVs using transpose() -> [ run_name, csv_file ]
        if (params.get('map_file', false) && params.get('truth_bed', false)) {
            eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
            
            // Combine each evaluation with its calling target BED file
            eval_joined = eval_inputs.combine(target_beds_ch, by: 0)
            
            EVALUATE_ACMG_CSV(
                eval_joined.map { [ it[0], it[1], it[2] ] }, // [ run_name, csv_file, bed_file ]
                file(params.map_file),
                file(params.truth_bed)
            )

            // 8. Dynamic Performance Comparison plotting
            PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
        }
    }
}
