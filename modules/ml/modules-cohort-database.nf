#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// COHORT FEATURE DATABASE MODULE
// Aggregates cohort-level consensus VCFs and sample-level extracted feature tables
// into high-performance relational and columnar databases (SQLite, DuckDB, Parquet)
// for machine learning model training and variant prioritization.
// =====================================================================================

process BUILD_COHORT_FEATURE_DATABASE {
    label 'feature_extraction'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort/database", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(cohort_vcf)
    path sample_features
    path ref_fasta
    path target_bed

    output:
    path "cohort_cnv_features.db", emit: sqlite_db
    path "cohort_cnv_features.variants.tsv.gz", emit: variants_tsv
    path "cohort_cnv_features.carriers.tsv.gz", emit: carriers_tsv
    path "cohort_cnv_features.parquet", emit: parquet_db, optional: true

    script:
    def fasta_arg = ref_fasta && ref_fasta.name != 'NO_FASTA' ? "--fasta ${ref_fasta}" : ""
    def bed_arg = target_bed && target_bed.name != 'NO_BED' ? "--target-bed ${target_bed}" : ""
    """
    python3 ${projectDir}/bin/build_cohort_feature_database.py \\
        --cohort-vcf "${cohort_vcf}" \\
        --sample-features-list ${sample_features} \\
        ${fasta_arg} \\
        ${bed_arg} \\
        --output-prefix "cohort_cnv_features"
    """
}
