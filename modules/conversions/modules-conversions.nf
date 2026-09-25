#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * MODULE: modules-conversions.nf
 *
 * Implements verified alignment format conversions between CRAM and BAM
 * using the 44-Point Forensic Verification Suite:
 *   - CONVERT_CRAM_TO_BAM_VERIFIED
 *   - CONVERT_BAM_TO_CRAM_VERIFIED
 */

process CONVERT_CRAM_TO_BAM_VERIFIED {
    tag "$sample_id"
    label 'samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/converted_bams", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(cram), path(crai)
    path fasta
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam_and_index

    script:
    """
    bash ${projectDir}/bin/convert_cram_to_bam_verified.sh \\
        ${cram} \\
        ${fasta} \\
        ${sample_id}.bam \\
        ${task.cpus}
    """
}

process CONVERT_CRAM_TO_TARGET_PADDED_BAM {
    tag "$sample_id"
    label 'samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/converted_bams", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(cram), path(crai)
    path fasta
    path fasta_fai
    path target_bed

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam_and_index

    script:
    def pad = params.get('padding_bp', 250)
    """
    bash ${projectDir}/bin/convert_cram_to_target_padded_bam.sh \\
        ${cram} \\
        ${fasta} \\
        ${target_bed} \\
        ${sample_id}.bam \\
        ${task.cpus} \\
        ${pad}
    """
}

process CONVERT_BAM_TO_CRAM_VERIFIED {
    tag "$sample_id"
    label 'samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/converted_crams", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}.cram"), path("${sample_id}.cram.crai"), emit: cram_and_index

    script:
    """
    bash ${projectDir}/bin/convert_bam_to_cram_verified.sh \\
        ${bam} \\
        ${fasta} \\
        ${sample_id}.cram \\
        ${task.cpus}
    """
}

workflow RUN_CRAM_TO_BAM_CONVERSION {
    take:
    cram_pairs_ch   // tuple val(sample_id), path(cram), path(crai)
    fasta
    fasta_fai

    main:
    CONVERT_CRAM_TO_BAM_VERIFIED(cram_pairs_ch, fasta, fasta_fai)

    emit:
    bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
}

workflow RUN_BAM_TO_CRAM_CONVERSION {
    take:
    bam_pairs_ch    // tuple val(sample_id), path(bam), path(bai)
    fasta
    fasta_fai

    main:
    CONVERT_BAM_TO_CRAM_VERIFIED(bam_pairs_ch, fasta, fasta_fai)

    emit:
    crams = CONVERT_BAM_TO_CRAM_VERIFIED.out.cram_and_index
}

workflow RUN_CRAM_TO_TARGET_PADDED_BAM_CONVERSION {
    take:
    cram_pairs_ch   // tuple val(sample_id), path(cram), path(crai)
    fasta
    fasta_fai
    target_bed

    main:
    CONVERT_CRAM_TO_TARGET_PADDED_BAM(cram_pairs_ch, fasta, fasta_fai, target_bed)

    emit:
    bams = CONVERT_CRAM_TO_TARGET_PADDED_BAM.out.bam_and_index
}
