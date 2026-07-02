#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// PROCESS: ANALYZE_CNV_QUALITY_METRICS (Avg depth, GC, Mappability, count per sample + Plots)
// =====================================================================================
process ANALYZE_CNV_QUALITY_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/metrics_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
    path probes_bed

    output:
    tuple val(run_name), path("${run_name}_cnv_metrics_summary.txt"), path("${run_name}_cnv_vcf_metrics.csv"), emit: metrics_reports
    path "qc_plots/cnv_qc_dashboard.png", emit: qc_dashboard

    script:
    """
    analyze_cnv_quality_metrics.py \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --vcf-dir . \\
        --probes-bed "${probes_bed}" \\
        --run-name "${run_name}" \\
        --output-summary "${run_name}_cnv_metrics_summary.txt" \\
        --output-csv "${run_name}_cnv_vcf_metrics.csv"

    plot_cnv_qc_metrics.py \\
        --metrics-csv "${run_name}_cnv_vcf_metrics.csv" \\
        --output-dir qc_plots
    """
}

// =====================================================================================
// SUB-WORKFLOW: CNV_QC_METRIC_EVALUATION
// =====================================================================================
workflow CNV_QC_METRIC_EVALUATION {
    take:
        annotated_ch // Channel: tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
        probes_bed   // Path: standard Sanger design target BED file

    main:
        ANALYZE_CNV_QUALITY_METRICS(annotated_ch, probes_bed)

    emit:
        metrics_reports = ANALYZE_CNV_QUALITY_METRICS.out.metrics_reports
        qc_dashboard    = ANALYZE_CNV_QUALITY_METRICS.out.qc_dashboard
}
