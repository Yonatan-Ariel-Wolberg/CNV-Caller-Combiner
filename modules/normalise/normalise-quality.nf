#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// PROCESS: NORMALISE_CNV_QUALITY_SCORES
// Normalise CNV quality scores to a common scale.
// Executed strictly downstream of quality score filtering (or on raw calls if
// quality filtering is skipped).
// =====================================================================================

process NORMALISE_CNV_QUALITY_SCORES {
    tag "${vcf.simpleName}"
    label 'pysam'
    publishDir { "${params.outdir}/${dir_suffix}/vcfs" }, mode: 'copy', overwrite: true

    input:
    path vcf
    val  caller_name  // value passed to --caller (e.g. 'CANOES', 'GATK')
    val  dir_suffix   // publish directory suffix  (e.g. 'out_CANOES', 'out_GCNV')

    output:
    path("*.normalised*.cnv.vcf.gz"),     emit: normalised_vcf
    path("*.normalised*.cnv.vcf.gz.tbi"), emit: normalised_vcf_index

    script:
    def caller = caller_name.toLowerCase()
    if (caller == 'gatk') caller = 'gcnv'
    def base_id = vcf.name
        .replaceAll(/\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth).*/, '')
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*/, '')
        .replaceAll(/(\.sorted|\.qs-filtered|\.filtered|\.cnv|\.annotated)*\.(vcf|vcf\.gz)$/, '')
        .replaceAll(/\.md$/, '')
    def is_filtered = vcf.name.contains('.qs-filtered') || vcf.name.contains('_filtered') || vcf.name.contains('.filtered')
    def normalised_base = is_filtered
        ? "${base_id}.${caller}.normalised.qs-filtered.cnv"
        : "${base_id}.${caller}.normalised.cnv"
    def normalised_gz   = "${normalised_base}.vcf.gz"
    """
    normalise_cnv_caller_quality_scores.py \\
        --input_vcf ${vcf} \\
        --output_vcf ${normalised_base}.vcf \\
        --caller ${caller_name}
    python3 -c "import pysam; pysam.tabix_index('${normalised_base}.vcf', preset='vcf', force=True)"
    rm -f ${normalised_base}.vcf
    """
}

// Sub-workflow wrapper for normalising a channel of VCFs
workflow NORMALISE {
    take:
        vcf_ch // channel of [ vcf, caller_name, dir_suffix ]

    main:
        NORMALISE_CNV_QUALITY_SCORES(vcf_ch)

    emit:
        normalised_vcf       = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf
        normalised_vcf_index = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf_index
}
