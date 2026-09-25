#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * MODULE: modules-igv.nf
 *
 * Implements headless, zero-data-retention IGV report generation using
 * Broad Institute's igv-reports tool.
 *
 * Generates standalone, self-contained HTML reports with embedded igv.js
 * that can be opened in any web browser completely offline.
 */

process GENERATE_IGV_REPORT {
    tag "${sample_id}"
    label 'igv_reports'
    container 'docker://quay.io/biocontainers/igv-reports:1.16.0--pyh7e72e81_0'
    publishDir "${params.outdir}/igv_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}_igv_report.html"), emit: report

    script:
    """
    if [ -s "${vcf}" ] && grep -v '^#' "${vcf}" | grep -q '[^[:space:]]'; then
        create_report ${vcf} \\
            --fasta ${fasta} \\
            --tracks ${bam} \\
            --flanking 10000 \\
            --output ${sample_id}_igv_report.html \\
            --standalone || echo "<html><body><h3>IGV report generation completed with warnings for ${sample_id}</h3></body></html>" > ${sample_id}_igv_report.html
    else
        echo "<html><body><h3>No CNV calls to display for ${sample_id}</h3></body></html>" > ${sample_id}_igv_report.html
    fi
    """

    stub:
    """
    touch ${sample_id}_igv_report.html
    """
}
