#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF } from '../common/bgzip-sort-index.nf'

// =====================================================================================
// CLAMMS MODULE
// =====================================================================================

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params.ref           = ''
params.probes        = ''
params.interval_list = ''
params.mappability   = ''
params.special_reg   = ''
params.sexinfo       = ''

// Define file variables for inputs
ref           = file(params.ref           ?: 'NO_FILE')
probes        = file(params.probes        ?: 'NO_FILE')
interval_list = file(params.interval_list ?: 'NO_FILE')
mappability   = file(params.mappability   ?: 'NO_FILE')
special_reg   = file(params.special_reg   ?: 'NO_FILE')
sexinfo       = file(params.sexinfo       ?: 'NO_FILE')
outdir        = file(params.outdir)

// =====================================================================================
// PROCESSES FOR CLAMMS
// =====================================================================================

// Process to generate windows for analysis
process GENERATE_WINDOWS {
    tag "Generate_Windows"
    label 'clamms_bedtools'

    output:
    path "windows.bed", emit: windows
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    export INSERT_SIZE=200
    sort -k1,1 -k2,2n ${probes} > targets_sorted.bed
    sort -k1,1 -k2,2n ${mappability} > mappability_sorted.bed
    \$CLAMMS_DIR/annotate_windows.sh targets_sorted.bed ${ref} mappability_sorted.bed \$INSERT_SIZE ${special_reg} > windows.bed
    """
}

// Process to calculate depth of coverage using samtools
process SAMTOOLS_DOC {
    tag "${sample_id}"
    label 'clamms_bedtools'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path windows

    output:
    tuple val(sample_id), path("${sample_id}.coverage.bed"), emit: coverage
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    samtools bedcov -Q 30 ${windows} ${bam} | awk '{ printf "%s\\t%d\\t%d\\t%.6g\\n", \$1, \$2, \$3, \$NF/(\$3-\$2); }' > ${sample_id}.coverage.bed
    """
}

// Process to normalize depth of coverage
process NORMALIZE_DOC {
    tag "${sample_id}"
    label 'clamms_bedtools'

    input:
    tuple val(sample_id), path(coverage)
    path windows

    output:
    tuple val(sample_id), path("${sample_id}.norm.cov.bed"), emit: norm_coverage
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    \$CLAMMS_DIR/normalize_coverage ${coverage} ${windows} | sed 's/^chr//g' > ${sample_id}.norm.cov.bed
    """
}

// Process to create PCA input data
process CREATE_PCA_DATA {
    tag "PCA_Data"
    label 'clamms_bedtools'

    input:
    path norm_covs

    output:
    path "pca.coordinates.txt", emit: pca_data
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    custom_ref_panel.sh
    """
}

// Process to get QC metrics using Picard
process GET_PICARD_QC_METRICS {
    tag "${sample_id}"
    label 'picard'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hs_metrics.txt"), emit: qc_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Collecting QC metrics for ${sample_id}..."
    picard CollectHsMetrics \\
        -I ${bam} \\
        -O ${sample_id}.hs_metrics.txt \\
        -R ${ref} \\
        -BI ${interval_list} \\
        -TI ${interval_list}
    echo "QC metrics collected: ${sample_id}.hs_metrics.txt"
    """
}

// Process to get mean insert size metrics using Picard
process GET_PICARD_MEAN_INSERT_SIZE {
    tag "${sample_id}"
    label 'picard'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.insert_size_metrics.txt"), emit: ins_size_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Collecting mean insert size metrics for ${sample_id}..."
    picard CollectInsertSizeMetrics \\
        -I ${bam} \\
        -O ${sample_id}.insert_size_metrics.txt \\
        -H ${sample_id}.insert_size_histogram.pdf
    echo "Insert size metrics collected: ${sample_id}.insert_size_metrics.txt"
    """
}

// Process to combine all QC metrics from Picard
process COMBINE_PICARD_QC_METRICS {
    tag "Combine_QC"
    label 'picard'

    input:
    path metrics

    output:
    path "qcs_metrics", emit: qcs_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    combine_picard_qc_metrics.sh
    """
}

// Process to create a custom reference panel
process CREATE_CUSTOM_REF_PANEL {
    tag "Ref_Panel"
    label 'R'

    input:
    path norm_covs
    path pca_data
    path combined_qc
    
    output:
    path "*.pdf", emit: plots
    path "*.ref.panel.files.txt", emit: ref_panel
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    custom_ref_panel.R ${pca_data} ${combined_qc} ${sexinfo}
    """
}

// Process to train models using the reference panel
process TRAIN_MODELS {
    tag "${sample_id}"
    label 'clamms_bedtools'

    input:
    tuple val(sample_id), path(ref_panel_file)
    path windows
    path norm_covs

    output:
    tuple val(sample_id), path("${sample_id}.models.bed"), emit: sample_models
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    sed 's/chr//g' ${windows} > windows.new.bed
    \$CLAMMS_DIR/fit_models ${ref_panel_file} windows.new.bed > ${sample_id}.models.bed
    """
}

// Process to call CNVs from the trained models
process CALL_CNVS {
    tag "${sample_id}"
    label 'clamms_bedtools'
    publishDir "${outdir}/out_CLAMMS/calls", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(sample_models), path(norm_cov)

    output:
    tuple val(sample_id), path("${sample_id}.cnv.bed"), emit: cnvs
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    sex=\$(grep "^${sample_id}\t" ${sexinfo} | cut -f 2)
    if [ -z "\$sex" ]; then echo "Error: Sample ${sample_id} not found in ${sexinfo}" >&2; exit 1; fi
    if ! \$CLAMMS_DIR/call_cnv ${norm_cov} ${sample_models} --sex \$sex > ${sample_id}.cnv.bed; then
        echo "WARNING: CLAMMS call_cnv failed for sample ${sample_id}. Creating empty BED file." >&2
        touch ${sample_id}.cnv.bed
    fi
    """
}

// Process to combine all sample CNV BED files into a unified cohort BED
// (Quality score filtering is handled centrally in quality-filters.nf)
process COMBINE_CLAMMS_CNVS {
    tag "Combine_CLAMMS"
    label 'clamms_bedtools'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    path cnvs

    output:
    path("samples.cnv.bed"), emit: all_cnvs
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    cat *.cnv.bed > samples.cnv.bed
    """
}

// Process to convert CLAMMS BED to VCF using the Python script
process CONVERT_CLAMMS_TO_VCF {
    tag "BED_TO_VCF"
    label 'clamms_bedtools'
    publishDir "${outdir}/out_CLAMMS/vcfs", mode: 'copy', overwrite: true

    input:
    path input_file
    path sample_file
    path fai_file

    output:
    path("*_CLAMMS_output.vcf"), emit: vcfs

    script:
    """
    clamms_bed_to_vcf.py --input_file ${input_file} --sample_file ${sample_file} --output_dir . --fai_file ${fai_file}
    echo "Conversion from BED to VCF completed"
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow CLAMMS {
    take:
    bam_ch         // channel: tuple(sample_id, bam, bai)
    fai_ch         // path: reference FASTA index (.fai)
    sample_file_ch // path: sample file for VCF conversion

    main:
    // Step 1: Generate analysis windows
    GENERATE_WINDOWS()

    // Step 2: Calculate depth of coverage per sample
    SAMTOOLS_DOC(bam_ch, GENERATE_WINDOWS.out.windows)

    // Step 3: Normalize depth of coverage per sample
    NORMALIZE_DOC(SAMTOOLS_DOC.out.coverage, GENERATE_WINDOWS.out.windows)

    // Step 4: Collect normalized coverages and create PCA data
    norm_coverage_ch = NORMALIZE_DOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }
        .ifEmpty { error "CLAMMS: no normalized coverage files were produced. Check upstream BAM/sample inputs." }
    CREATE_PCA_DATA(norm_coverage_ch.collect())

    // Step 5: Get Picard QC metrics per sample
    GET_PICARD_QC_METRICS(bam_ch)
    GET_PICARD_MEAN_INSERT_SIZE(bam_ch)

    // Step 6: Combine all Picard QC metrics
    all_metrics_ch = GET_PICARD_QC_METRICS.out.qc_metrics
        .mix(GET_PICARD_MEAN_INSERT_SIZE.out.ins_size_metrics)
        .map { sample_id, metrics -> metrics }
        .collect()

    COMBINE_PICARD_QC_METRICS(all_metrics_ch)

    // Step 7: Create custom reference panel using all normalized coverages
    CREATE_CUSTOM_REF_PANEL(
        NORMALIZE_DOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }.collect(),
        CREATE_PCA_DATA.out.pca_data,
        COMBINE_PICARD_QC_METRICS.out.qcs_metrics
    )

    // Step 8: Train models per sample using its reference panel file
    ref_panel_per_sample_ch = CREATE_CUSTOM_REF_PANEL.out.ref_panel
        .flatten()
        .map { file -> tuple(file.name.replace('.ref.panel.files.txt', ''), file) }

    TRAIN_MODELS(
        ref_panel_per_sample_ch,
        GENERATE_WINDOWS.out.windows,
        NORMALIZE_DOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }.collect()
    )

    // Step 9: Call CNVs per sample by joining models with normalized coverage
    call_cnv_input_ch = TRAIN_MODELS.out.sample_models
        .join(NORMALIZE_DOC.out.norm_coverage, by: 0)

    CALL_CNVS(call_cnv_input_ch)

    // Step 10: Collect all per-sample CNV BED calls into unified cohort BED
    // (Quality score filtering is handled centrally in quality-filters.nf)
    clamms_cnvs_ch = CALL_CNVS.out.cnvs.map { sample_id, cnv_file -> cnv_file }
        .ifEmpty { error "CLAMMS: no CNV calls were produced. Check workflow inputs and model training outputs." }
    COMBINE_CLAMMS_CNVS(clamms_cnvs_ch.collect())

    // Step 11: Convert raw combined BED to VCF
    CONVERT_CLAMMS_TO_VCF(
        COMBINE_CLAMMS_CNVS.out.all_cnvs,
        sample_file_ch,
        fai_ch
    )

    // Step 12: Sort, compress, and index the VCF
    BGZIP_SORT_INDEX_VCF(CONVERT_CLAMMS_TO_VCF.out.vcfs.flatten(), 'CLAMMS', 'out_CLAMMS')

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
