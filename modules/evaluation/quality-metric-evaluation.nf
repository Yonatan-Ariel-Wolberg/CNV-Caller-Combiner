#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================================
//  QUALITY METRIC EVALUATION MODULE
//  Extracts, assesses, and visualizes CNV callset quality metrics and variant
//  attrition dynamics across all pipeline stages and strata:
//    1. Variant counts, DEL / DUP counts & DEL/DUP ratio
//    2. Length distributions (Median bp, IQR, Total Genomic Footprint in Mb)
//    3. Target probe support (Median probes/call, single-probe vs multi-probe %)
//    4. Quality score distributions (QUAL)
//    5. Sample-level distribution (Mean/median calls per sample, outlier detection)
//    6. Target depth of coverage, GC content, and mappability analysis
//    7. Multi-stage summary reports (.md, .tsv)
//    8. Variant Attrition Funnel tracking candidate progression & drop-off
// ================================================================================

process ANALYZE_CNV_QUALITY_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_metrics", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
    path probes_bed

    output:
    tuple val(run_name), path("${run_name}_cnv_metrics_summary.txt"), path("${run_name}_cnv_vcf_metrics.csv"), emit: metrics_reports

    script:
    """
    python3 "${projectDir}/bin/analyze_cnv_quality_metrics.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --vcf-dir . \\
        --probes-bed "${probes_bed}" \\
        --run-name "${run_name}" \\
        --output-summary "${run_name}_cnv_metrics_summary.txt" \\
        --output-csv "${run_name}_cnv_vcf_metrics.csv"
    """
}

process GENERATE_ATTRITION_FUNNEL {
    label 'feature_extraction'
    publishDir "${params.outdir}/attrition_funnel", mode: 'copy', overwrite: true

    input:
    path eval_reports

    output:
    path "attrition_funnel_table.md", emit: table
    path "attrition_funnel.png", emit: plot

    script:
    """
    python3 "${projectDir}/bin/generate_attrition_funnel.py" \\
        --results-dir "${params.outdir}" \\
        --output-prefix attrition_funnel
    """
}

process ASSESS_CALLSET_QUALITY_METRICS {
    tag "${label}"
    label 'feature_extraction'
    publishDir "${params.outdir}/callset_quality_assessments/${label}", mode: 'copy', overwrite: true

    input:
    tuple val(label), path(vcfs)
    path probes_bed

    output:
    tuple val(label), path("${label}_quality_summary.txt"), path("${label}_quality_metrics.csv"), emit: report
    path "${label}_quality_summary.txt", emit: summary_file

    script:
    def probes_arg = (probes_bed && probes_bed.name != 'NO_BED' && probes_bed.exists()) ? "--probes-bed ${probes_bed}" : ""
    """
    python3 "${projectDir}/bin/assess_callset_quality_metrics.py" \\
        --vcf-files ${vcfs} \\
        --stage-name "${label}" \\
        ${probes_arg} \\
        --output-summary "${label}_quality_summary.txt" \\
        --output-csv "${label}_quality_metrics.csv"
    """
}

process COMPILE_ALL_STAGES_QUALITY_REPORT {
    label 'feature_extraction'
    publishDir "${params.outdir}/callset_quality_assessments", mode: 'copy', overwrite: true

    input:
    path summary_reports

    output:
    path "all_stages_quality_metrics_summary.md", emit: markdown_table
    path "all_stages_quality_metrics_summary.tsv", emit: tsv_table

    script:
    """
    python3 "${projectDir}/bin/compile_quality_summary.py" \\
        --output-md all_stages_quality_metrics_summary.md \\
        --output-tsv all_stages_quality_metrics_summary.tsv
    """
}
