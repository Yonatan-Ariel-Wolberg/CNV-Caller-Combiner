#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESS: NORMALISE CNV QUALITY SCORES (standalone, caller-agnostic)
// =====================================================================================
//
// Runs bin/normalise_cnv_caller_quality_scores.py on a pre-existing caller VCF.
//
// This process is the standalone equivalent of the NORMALISE_CNV_QUALITY_SCORES
// process that is embedded inside each caller module (CANOES, CLAMMS, XHMM,
// GATK-gCNV, CNVkit, DRAGEN, INDELIBLE).  It allows the normalisation step –
// and therefore the qual_norm_{caller} features used by the feature-extraction
// workflow – to be produced even when the full caller pipeline has not been run
// (e.g. when pre-existing caller VCFs are supplied from an external source).
//
// Input tuple layout:
//   val  sample_name  -- bare sample identifier (no extension)
//   path vcf          -- raw caller VCF (plain .vcf or bgzipped .vcf.gz)
//   val  caller       -- one of: CANOES CLAMMS XHMM GATK CNVKIT DRAGEN INDELIBLE
//
// Outputs:
//   normalised_vcf       -- bgzipped normalised VCF  (*.normalised.vcf.gz)
//   normalised_vcf_index -- tabix index              (*.normalised.vcf.gz.tbi)

process NORMALISE_CNV_QUALITY_SCORES {
    tag "${sample_name} [${caller}]"
    label 'pysam'
    publishDir "${outdir}/out_NORMALISED/${caller}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_name), path(vcf), val(caller)

    output:
    tuple val(sample_name), path("*.normalised.vcf.gz"),     emit: normalised_vcf
    tuple val(sample_name), path("*.normalised.vcf.gz.tbi"), emit: normalised_vcf_index

    script:
    // Strip .sorted.vcf.gz / .vcf.gz / .vcf to get a clean base name for output files
    def base     = vcf.name.replaceAll(/(?i)\.(sorted\.)?vcf(\.gz)?$/, '')
    def norm_vcf = "${base}.normalised.vcf"
    def norm_gz  = "${norm_vcf}.gz"
    """
    normalise_cnv_caller_quality_scores.py \\
        --input_vcf  '${vcf}'      \\
        --output_vcf '${norm_vcf}' \\
        --caller     '${caller}'
    bgzip -c '${norm_vcf}' > '${norm_gz}'
    tabix -p vcf '${norm_gz}'
    """
}

// =====================================================================================
// SUB-WORKFLOW
// =====================================================================================

workflow NORMALISE {
    take:
        // Channel of tuples: [ sample_name, vcf_path, caller_name ]
        // caller_name must be one of the values accepted by
        // normalise_cnv_caller_quality_scores.py --caller:
        //   CANOES  CLAMMS  XHMM  GATK  CNVKIT  DRAGEN  INDELIBLE
        vcf_ch

    main:
        NORMALISE_CNV_QUALITY_SCORES(vcf_ch)

    emit:
        normalised_vcf       = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf
        normalised_vcf_index = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf_index
}
