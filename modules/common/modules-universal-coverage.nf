#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * =====================================================================================
 * UNIVERSAL MOSDEPTH COVERAGE MODULE
 * =====================================================================================
 * Extracts per-target coverage ONCE per sample using Mosdepth and converts the coverage
 * directly into caller-ready inputs for:
 *   - XHMM:     DATA.RD.txt matrix
 *   - CANOES:   Per-chromosome ${chr}_canoes_reads.txt and sample lists
 *   - CLAMMS:   ${sample}.coverage.bed
 *   - CNVkit:   ${sample}.targetcoverage.cnn
 *
 * Eliminates redundant BAM parsing by 4 independent callers, cutting ~80% of cluster I/O.
 * =====================================================================================
 */

process EXTRACT_MOSDEPTH_COVERAGE {
    tag "$sample_id"
    label 'mosdepth'
    cpus 2
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fasta_fai
    path target_bed

    output:
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"),     emit: regions
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary, optional: true

    script:
    def is_cram = bam.name.endsWith('.cram')
    def fasta_arg = is_cram ? "--fasta ${fasta}" : ""
    """
    mosdepth -t ${task.cpus} \\
             -n \\
             -Q 20 \\
             --by ${target_bed} \\
             ${fasta_arg} \\
             ${sample_id} \\
             ${bam}
    """
}

process CONVERT_MOSDEPTH_TO_XHMM {
    tag "xhmm_matrix"
    label 'cnvkit'
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path regions_beds
    val sample_ids

    output:
    path "DATA.RD.txt", emit: xhmm_matrix

    script:
    def samples_arg = sample_ids instanceof List ? sample_ids.join(' ') : sample_ids
    """
    python3 ${projectDir}/bin/convert_mosdepth_to_caller_inputs.py xhmm \\
        --beds ${regions_beds} \\
        --samples ${samples_arg} \\
        -o DATA.RD.txt
    """
}

process CONVERT_MOSDEPTH_TO_CANOES {
    tag "canoes_matrices"
    label 'cnvkit'
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}/out_CANOES", mode: 'copy', overwrite: true

    input:
    path regions_beds
    val sample_ids

    output:
    path "canoes_counts/*_canoes_reads.txt", emit: count_matrices
    path "canoes_counts/*_sample_list",      emit: sample_lists

    script:
    def samples_arg = sample_ids instanceof List ? sample_ids.join(' ') : sample_ids
    """
    python3 ${projectDir}/bin/convert_mosdepth_to_caller_inputs.py canoes \\
        --beds ${regions_beds} \\
        --samples ${samples_arg} \\
        -o canoes_counts
    """
}

process CONVERT_MOSDEPTH_TO_CLAMMS {
    tag "$sample_id"
    label 'cnvkit'
    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id), path(regions_bed)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.bed"), emit: coverage

    script:
    """
    python3 ${projectDir}/bin/convert_mosdepth_to_caller_inputs.py clamms \\
        --bed ${regions_bed} \\
        -o ${sample_id}.coverage.bed
    """
}

process CONVERT_MOSDEPTH_TO_CNVKIT {
    tag "$sample_id"
    label 'cnvkit'
    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id), path(regions_bed)
    path target_bed

    output:
    tuple val(sample_id), path("${sample_id}.targetcoverage.cnn"), emit: target_cov

    script:
    """
    python3 ${projectDir}/bin/convert_mosdepth_to_caller_inputs.py cnvkit \\
        --bed ${regions_bed} \\
        --target-bed ${target_bed} \\
        -o ${sample_id}.targetcoverage.cnn
    """
}
