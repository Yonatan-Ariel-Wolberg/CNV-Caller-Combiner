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
params.precalled_control_vcfs  = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/control/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-3264cf37-6594-4c0c-b283-6beecac654cf/EGAF*/*.cnv.vcf.gz"
params.precalled_mosdepth_vcfs = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/consensus/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-8dc80f47-9ded-4d22-9d9b-938dfb5bcabd/EGAF*/*.cnv.vcf.gz"
params.precalculated_consensus_bed = "/home/ywolberg/dragen-germline-enrichment-v4-4-6/multiqc-filter-mosdepth/ddduk_karyotypic_consensus.bed"

// 5. Cohort Filtering & ACMG Annotations
params.total_samples         = 90
params.max_af                = 0.05
params.outdir                = "./results"
params.annotsv_db            = ""
params.annotsv_version       = "3.5.10"

// 6. Evaluator Paths
params.map_file              = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/DDDUK-MAP/EGAF-EGAN-ddduk.map"
params.truth_bed             = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed"

// Define dynamic repository and bin directory paths depending on script location
def repoDir = file("${baseDir}/main.nf").exists() ? "${baseDir}" : (file("${baseDir}/../main.nf").exists() ? "${baseDir}/.." : "${baseDir}")
def binDir = file("${repoDir}/bin").exists() ? "${repoDir}/bin" : "${repoDir}"


/*
 * =====================================================================================
 *  PHASE 1: MOSDEPTH TARGET BED GENERATION
 * =====================================================================================
 */

process RUN_MOSDEPTH {
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0'
    cpus 4
    memory '24 GB'
    maxForks 10
    
    input:
    tuple val(sample), path(cram), path(crai)
    path fasta
    path fasta_fai
    path bed

    output:
    path "${sample}.regions.bed.gz", emit: usable_beds
    
    script:
    """
    mosdepth --threads ${task.cpus} \\
             --fasta ${fasta} \\
             --by ${bed} \\
             ${sample} ${cram}
    """
}

process BUILD_CONSENSUS {
    container 'quay.io/biocontainers/bedtools:2.27.1--h077b44d_9'
    cpus 2
    memory '8 GB'
    
    publishDir "${params.outdir}/consensus_bed", mode: 'copy'

    input:
    path usable_beds
    path bed

    output:
    path "ddduk_karyotypic_consensus.bed", emit: karyotypic_bed

    script:
    """
    TOTAL_SAMPLES=\$(ls *.regions.bed.gz | wc -l)
    MIN_PASS=\$(( (TOTAL_SAMPLES * ${params.pass_percent}) / 100 ))

    # Extract target interval coverages for each sample
    for f in *.regions.bed.gz; do
        gunzip -c "\$f" | cut -f4 > tmp_\$f
    done
    
    # Get coordinates from the first file
    first_file=\$(ls *.regions.bed.gz | head -n1)
    gunzip -c "\$first_file" | cut -f1-3 > coords.txt
    
    # Combine coordinates and coverages
    paste coords.txt tmp_*.regions.bed.gz > matrix.txt
    
    # Filter intervals based on sample coverage count
    awk -v min="\$MIN_PASS" -v n="\$TOTAL_SAMPLES" 'BEGIN{FS="\\t"; OFS="\\t"} {
        count = 0;
        for(i=4; i<=3+n; i++) {
            if(\$i >= 20.0) count++;
        }
        if(count >= min) {
            print \$1, \$2, \$3;
        }
    }' matrix.txt > ddduk_clean_targeted_consensus.bed
    
    # Clean up temp files and sort consensus bed
    rm -f tmp_* matrix.txt coords.txt
    sort -k1,1V -k2,2n ddduk_clean_targeted_consensus.bed | bedtools merge -i stdin -d 1500 > ddduk_karyotypic_consensus.bed
    """
}


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

process bcftoolsFilterVcfs {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*_filtered.cnv.vcf'), emit: filtered_data

    script:
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.cnv.vcf.gz; do
        base=\$(basename "\$vcf" .cnv.vcf.gz)
        filtered_vcf="\${base}_filtered.cnv.vcf"
        
        # Filter for quality PASS variants
        bcftools view -i 'FILTER="PASS" && QUAL>=15 && FORMAT/BC>=3' "\$vcf" \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                svtype="CNV";
                if (\$5 == "<DEL>") svtype="DEL";
                else if (\$5 == "<DUP>") svtype="DUP";
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                print;
            }' > "\$filtered_vcf"
            
        if grep -q -m1 -v '^#' "\$filtered_vcf"; then
            echo "\$filtered_vcf" >> vcf_list.txt
        else
            rm -f "\$filtered_vcf"
        fi
    done
    """
}

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

process bcftoolsFlattenCohort {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf)

    output:
    tuple val(run_name), path('cohort_cnvs_by_sample.tsv'), emit: cohort_tsv

    script:
    """
    bcftools query -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' cohort_merged.vcf \\
        | awk '\$6 != "0/0" && \$6 != "./." && \$6 != "0"' \\
        > cohort_cnvs_by_sample.tsv
    """
}

process filterRareCnvs {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(cohort_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), path('rare_regions.tsv'), emit: rare_data

    script:
    """
    # cache invalidate for small cohort variant filtering fix
    python3 "${binDir}/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --total-samples ${params.total_samples} \\
        --max-af ${params.max_af} \\
        --output-regions "rare_regions.tsv"
    """
}

process filterVcfByRareRegions {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf), path(rare_regions)

    output:
    tuple val(run_name), path('cohort_rare.vcf'), emit: rare_vcf

    script:
    """
    if [ ! -s "${rare_regions}" ]; then
        bcftools view -h "${merged_vcf}" -o cohort_rare.vcf
    else
        bcftools view -T "${rare_regions}" "${merged_vcf}" -o cohort_rare.vcf
    fi
    """
}

process DOWNLOAD_ANNOTSV_DB {
    publishDir "${repoDir}/annotsv_db", mode: 'copy'

    output:
    path "AnnotSV_annotations", emit: db_dir

    script:
    def version_tag = params.annotsv_version.startsWith('v') ? params.annotsv_version : "v${params.annotsv_version}"
    def version_clean = params.annotsv_version.replaceAll(/^v/, '')
    def parts = version_clean.split('\\.')
    def major_minor = parts.size() >= 2 ? "${parts[0]}.${parts[1]}" : version_clean
    """
    git clone --branch ${version_tag} --depth 1 https://github.com/lgmgeo/AnnotSV.git
    cd AnnotSV
    
    # Pre-download the correct human annotations file. Keep it under the name the Makefile expects
    # (Annotations_Human_${major_minor}.tar.gz) and copy it to the versioned clean name just in case.
    curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_${major_minor}.tar.gz
    cp Annotations_Human_${major_minor}.tar.gz Annotations_Human_${version_clean}.tar.gz 2>/dev/null || true
    
    make PREFIX=\$(pwd)/../AnnotSV_annotations install
    make PREFIX=\$(pwd)/../AnnotSV_annotations install-human-annotation
    cd ..
    rm -rf AnnotSV
    """
}

process runAnnotSV {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(rare_vcf)
    path db_dir

    output:
    tuple val(run_name), path('cohort_rare.annotated.tsv'), emit: annotated_tsv

    script:
    """
    export ANNOTSV=\$(readlink -f ${db_dir})
    
    # Run AnnotSV
    AnnotSV \\
        -SVinputFile "${rare_vcf}" \\
        -outputDir . \\
        -outputFile cohort_rare.annotated \\
        -genomeBuild GRCh38 \\
        -overwrite 0 || true
        
    # If the output file was not created (e.g., empty input), touch it with standard headers
    if [ ! -f cohort_rare.annotated.tsv ]; then
        echo -e "Annotation_mode\\tSV_chrom\\tSV_start\\tSV_end\\tSV_type\\tACMG_class\\tAnnotSV_ranking_score\\tAnnotSV_ranking_criteria\\tGene_name" > cohort_rare.annotated.tsv
    fi
    """
}

process addACMGAnnotations {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/acmg_reports", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), emit: rare_cnvs_annotated
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs_*.csv'), emit: score_separated_csvs

    script:
    """
    python3 "${binDir}/add-acmg-annotations.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv"
    """
}

process ANALYZE_CNV_QUALITY_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/metrics_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
    path probes_bed

    output:
    tuple val(run_name), path("${run_name}_cnv_metrics_summary.txt"), path("${run_name}_cnv_vcf_metrics.csv"), emit: metrics_reports

    script:
    """
    python3 "${binDir}/analyze_cnv_quality_metrics.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --vcf-dir . \\
        --probes-bed "${probes_bed}" \\
        --run-name "${run_name}" \\
        --output-summary "${run_name}_cnv_metrics_summary.txt" \\
        --output-csv "${run_name}_cnv_vcf_metrics.csv"
    """
}


/*
 * =====================================================================================
 *  PHASE 4: RIGOROUS EVALUATION & VISUALIZATION
 * =====================================================================================
 */

process EVALUATE_ACMG_CSV {
    tag "${run_name}-${csv.baseName}"
    publishDir "${params.outdir}/evaluate-results/${run_name}", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(csv), path(probes_bed)
    path map_file
    path truth_bed

    output:
    path "*_performance_report.txt", emit: report

    script:
    def base_name = csv.baseName
    """
    # cache invalidate for evaluation sample ID mapping fix
    python3 "${binDir}/evaluate_separated_helper.py" \\
        --input-csv "${csv}" \\
        --map-file "${map_file}" \\
        --truth-bed "${truth_bed}" \\
        --probes-bed "${probes_bed}" \\
        --workflow-type "${run_name}" \\
        --output "${run_name}_${base_name}_performance_report.txt"
    """
}

process PLOT_EVALUATION_METRICS {
    publishDir "${params.outdir}/evaluate-results", mode: 'copy', overwrite: true

    input:
    path reports
    
    output:
    path 'performance_metrics_comparison.png', emit: comparison_plot

    script:
    """
    # Structure files into subdirectories so parser detects Control vs Mosdepth
    mkdir -p Control Mosdepth
    for f in *_performance_report.txt; do
        if [[ "\$f" == *"Control"* ]]; then
            cp "\$f" Control/
        elif [[ "\$f" == *"Mosdepth"* ]]; then
            cp "\$f" Mosdepth/
        else
            cp "\$f" Control/
        fi
    done
    
    python3 "${binDir}/plot_evaluation_performance_dynamic.py" \\
        --input-dir . \\
        --output-png performance_metrics_comparison.png
    """
}


/*
 * =====================================================================================
 *  WORKFLOW ENTRY POINT
 * =====================================================================================
 */

params.annotsv_version       = "3.5.10"
params.skip_annotation       = false

// 6. Evaluator Paths
params.map_file              = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/DDDUK-MAP/EGAF-EGAN-ddduk.map"
params.truth_bed             = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed"

// Define dynamic repository and bin directory paths depending on script location
def repoDir = file("${baseDir}/main.nf").exists() ? "${baseDir}" : (file("${baseDir}/../main.nf").exists() ? "${baseDir}/.." : "${baseDir}")
def binDir = file("${repoDir}/bin").exists() ? "${repoDir}/bin" : "${repoDir}"


workflow {
    // 1. Reference Validation & Inputs Setup
    def modes = params.run_modes.tokenize(',')
    
    // Validate that inputs for mosdepth are present if Mosdepth run mode is selected
    if (modes.contains('Mosdepth')) {
        if (!params.bam_dir) error "DRAGEN consensus generation (Mosdepth mode) requires --bam_dir containing BAM files."
    }

    // Initialize the AnnotSV database path channel
    if (params.annotsv_db) {
        annotsv_db_ch = Channel.fromPath(params.annotsv_db)
    } else {
        DOWNLOAD_ANNOTSV_DB()
        annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
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
            if (file("${params.outdir}/${run_name}/raw_calls/out_vcfs/EGAF00001051389.md_filtered.cnv.vcf.gz").exists()) {
                log.info "Bypassing ICA DRAGEN calling for mode ${run_name} because raw calls were found in results directory."
                return false
            }
            return true
        }

        bypassed_vcfs_ch = target_beds_ch.filter { run_name, bed ->
            if (file("${params.outdir}/${run_name}/raw_calls/out_vcfs/EGAF00001051389.md_filtered.cnv.vcf.gz").exists()) {
                return true
            }
            return false
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
    bcftoolsFilterVcfs(vcfs_ch)

    if (!params.skip_annotation) {
        SURVIVORMergeCohort(bcftoolsFilterVcfs.out.filtered_data)
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
            .join(bcftoolsFilterVcfs.out.filtered_data.map { [ it[0], it[2] ] })
            
        // Combine with target BED files to get the correct bed file per run_name
        metrics_joined = metrics_inputs.combine(target_beds_ch, by: 0)
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_joined.map { [ it[0], it[1], it[2], it[3] ] }, // [run_name, rare_cnvs_csv, annotated_tsv, filtered_vcfs]
            metrics_joined.map { it[4] } // bed_file
        )

        // 7. Score-separated Evaluation
        // Expand the list of ACMG CSVs using transpose() -> [ run_name, csv_file ]
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
