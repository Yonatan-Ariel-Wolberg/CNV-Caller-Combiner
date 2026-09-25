#!/usr/bin/env nextflow

// =============================================================================
// Detailed documentation: clamms_germline_cnv_end_to_end.nf
// Purpose: This script supports CNV calling and/or post-processing tasks in
//   this project, with parameters and I/O conventions aligned to repository workflows.
// Scope: End-to-end run path for CLAMMS.
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
 *  CLAMMS GERMLINE CNV END-TO-END ANALYSIS PIPELINE
 *  Consensus BED Generation -> CLAMMS Caller -> Cohort Filtering & ACMG -> Evaluation & Plotting
 * -------------------------------------------------------------------------------------
 */

// 1. Reference & Input Paths
params.fasta                 = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.fasta_fai             = "${params.fasta}.fai"
params.dict                  = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.dict"
params.standard_bed          = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
params.control_interval_list = "/dataG/ddd/data/resources/canoes/probes_sanger.interval_list"
params.bam_dir               = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/bams"

// 2. Pipeline Modes & Options
params.run_modes             = "Control,Mosdepth" // Modes to run. Can be "Control", "Mosdepth", or "Control,Mosdepth"
params.run_calling           = true  // Execute CNV calling locally using CLAMMS
params.pass_percent          = 90    // Consensus threshold: regions covered in >= X% of samples

// 3. Pre-called VCF Fallback (used if params.run_calling = false)
params.precalled_control_vcfs  = "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/10-samples/clamms/out_CLAMMS/vcfs/*.normalised.vcf.gz"
params.precalled_mosdepth_vcfs = "" // path to pre-called mosdepth CLAMMS VCFs if run_calling = false
params.precalculated_consensus_bed = null

// 4. Cohort Filtering & ACMG Annotations
params.total_samples         = 10
params.max_af                = 0.05
params.outdir                = "./results_clamms_e2e"
params.annotsv_db            = ""

// 5. Evaluator Paths
params.map_file              = "sample_id.map"
params.truth_bed             = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed"

// Define dynamic repository and bin directory paths depending on script location
def repoDir = file("${baseDir}/main.nf").exists() ? "${baseDir}" : (file("${baseDir}/../main.nf").exists() ? "${baseDir}/.." : "${baseDir}")
def binDir = file("${repoDir}/bin").exists() ? "${repoDir}/bin" : "${repoDir}"


/*
 * =====================================================================================
 *  PHASE 1: MOSDEPTH TARGET BED GENERATION & INTERVAL LIST CONVERSION
 * =====================================================================================
 */





process CONVERT_BED_TO_INTERVAL_LIST {
    container 'broadinstitute/gatk:4.6.2.0'
    label 'gatk'
    
    input:
    path bed
    path dict

    output:
    path "ddduk_karyotypic_consensus.interval_list", emit: interval_list

    script:
    """
    gatk BedToIntervalList \\
         -I ${bed} \\
         -O ddduk_karyotypic_consensus.interval_list \\
         -SD ${dict}
    """
}


/*
 * =====================================================================================
 *  PHASE 2: CLAMMS NESTED CALLING WORKFLOW
 * =====================================================================================
 */

process RUN_CLAMMS_CALLING {
    tag "${run_name}"
    maxForks 1
    cpus { [4, (Runtime.runtime.availableProcessors() ?: 2)].min() }
    memory '16 GB'
    publishDir "${params.outdir}/${run_name}/raw_calls", mode: 'copy'

    input:
    tuple val(run_name), path(bed_file), path(interval_list)

    output:
    tuple val(run_name), path("out_vcfs/*.vcf.gz"), emit: vcfs

    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Resolve absolute paths of input files from task directory and copy to persistent folder
    mkdir -p ${repoDir}/.nested_clamms_work/${run_name}/inputs
    cp -L ${bed_file} ${repoDir}/.nested_clamms_work/${run_name}/inputs/targets.bed
    cp -L ${interval_list} ${repoDir}/.nested_clamms_work/${run_name}/inputs/targets.interval_list
    
    # Create and move to persistent workspace directory
    cd ${repoDir}/.nested_clamms_work/${run_name}
    
    # Generate the samplesheet dynamically listing all BAM files
    echo -e "SampleID\\tBAM\\tBAI" > samplesheet.tsv
    for bam in ${params.bam_dir}/*.bam; do
        sample_id=\$(basename "\$bam" .bam)
        sample_id_clean=\$(echo "\$sample_id" | sed 's/\\.md//g')
        if [ -f "\${bam}.bai" ]; then
            bai="\${bam}.bai"
        elif [ -f "\${bam%.bam}.bai" ]; then
            bai="\${bam%.bam}.bai"
        else
            bai=""
        fi
        echo -e "\$sample_id_clean\\t\$bam\\t\$bai" >> samplesheet.tsv
    done
    
    # Execute nested nextflow run for CLAMMS
    nextflow run ${repoDir}/main.nf \\
        -profile wits,small \\
        --workflow clamms \\
        -params-file ${repoDir}/params/ddd-uk/params-clamms-wits-ddd-uk.json \\
        --probes "${repoDir}/.nested_clamms_work/${run_name}/inputs/targets.bed" \\
        --interval_list "${repoDir}/.nested_clamms_work/${run_name}/inputs/targets.interval_list" \\
        --outdir ./tmp_out \\
        --samplesheet_bams "\$(pwd)/samplesheet.tsv" \\
        -resume
        
    # Move back to parent task directory and copy raw calls to output folder
    cd \$OLDPWD
    mkdir -p out_vcfs
    cp ${repoDir}/.nested_clamms_work/${run_name}/tmp_out/out_CLAMMS/vcfs/*.normalised.vcf.gz out_vcfs/
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

    // 1. Parse run modes to execute
    def modes = params.run_modes.split(',').collect { it.trim() }

    // Set up AnnotSV database channel (uses host pre-installed DB or downloads dynamically)
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
    
    // 2. Generate consensus BED dynamically if 'Mosdepth' run mode is active
    if (modes.contains('Mosdepth')) {
        if (params.precalculated_consensus_bed && file(params.precalculated_consensus_bed).exists()) {
            mosdepth_bed_ch = Channel.value(file(params.precalculated_consensus_bed))
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
        CONVERT_BED_TO_INTERVAL_LIST(mosdepth_bed_ch, file(params.dict))
        mosdepth_interval_list_ch = CONVERT_BED_TO_INTERVAL_LIST.out.interval_list
    } else {
        mosdepth_bed_ch = Channel.empty()
        mosdepth_interval_list_ch = Channel.empty()
    }

    // 3. Construct target_beds_ch pairing run_name with target BED and GATK interval_list
    control_bed_ch = modes.contains('Control') ? Channel.of([ 'Control', file(params.standard_bed), file(params.control_interval_list) ]) : Channel.empty()
    mosdepth_bed_entry_ch = modes.contains('Mosdepth') ? mosdepth_bed_ch.combine(mosdepth_interval_list_ch).map { bed, il -> [ 'Mosdepth', bed, il ] } : Channel.empty()
    
    target_beds_ch = control_bed_ch.concat(mosdepth_bed_entry_ch)

    // 4. CNV Variant Calling Phase (Dynamic Calling vs Precalled Fallback)
    if (params.run_calling) {
        called_modes_ch = target_beds_ch.filter { run_name, bed, il ->
            def existing = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            if (existing && existing.size() > 0) {
                log.info "Bypassing CLAMMS CNV calling for mode ${run_name} because raw calls were found in results directory."
                return false
            }
            return true
        }

        bypassed_vcfs_ch = target_beds_ch.filter { run_name, bed, il ->
            def existing = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            return (existing && existing.size() > 0)
        }.map { run_name, bed, il ->
            def files = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            return [ run_name, files ]
        }

        RUN_CLAMMS_CALLING(called_modes_ch)
        vcfs_ch = RUN_CLAMMS_CALLING.out.vcfs.concat(bypassed_vcfs_ch)
    } else {
        vcfs_ch = Channel.fromList(modes)
            .map { run_name ->
                def glob_path = (run_name == 'Control') ? params.precalled_control_vcfs : params.precalled_mosdepth_vcfs
                if (!glob_path) error "Pre-called VCF path is not defined for mode ${run_name}."
                def files = file(glob_path)
                return [ run_name, files ]
            }
    }

    // 5. Cohort Analysis & ACMG Priorities
    FILTER_CLAMMS_QUALITY(vcfs_ch)

    if (!params.skip_annotation) {
        SURVIVORMergeCohort(FILTER_CLAMMS_QUALITY.out.filtered_data)
        bcftoolsFlattenCohort(SURVIVORMergeCohort.out.merged_vcf)
        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
        
        slice_joined = SURVIVORMergeCohort.out.merged_vcf.join(filterRareCnvs.out.rare_data)
        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
        
        runAnnotSV(filterVcfByRareRegions.out.rare_vcf, annotsv_db_ch)
        
        acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })

        // 6. CNV Quality Metrics Assessment
        metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
            .join(runAnnotSV.out.annotated_tsv)
            .join(FILTER_CLAMMS_QUALITY.out.filtered_data.map { [ it[0], it[2] ] })
            
        metrics_joined = metrics_inputs.combine(target_beds_ch.map { [ it[0], it[1] ] }, by: 0)
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_joined.map { [ it[0], it[1], it[2], it[3] ] },
            metrics_joined.map { it[4] }
        )

        // 7. Score-separated Evaluation
        if (params.get('map_file', false) && params.get('truth_bed', false)) {
            eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
            eval_joined = eval_inputs.combine(target_beds_ch.map { [ it[0], it[1] ] }, by: 0)
            
            EVALUATE_ACMG_CSV(
                eval_joined.map { [ it[0], it[1], it[2] ] },
                file(params.map_file),
                file(params.truth_bed)
            )

            // 8. Dynamic Comparison plotting
            PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
        }
    }
}
