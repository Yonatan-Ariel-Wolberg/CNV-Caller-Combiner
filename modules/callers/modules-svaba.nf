#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_SV_TUPLE_VCF } from '../common/bgzip-sort-index.nf'

/**
 * SVABA Structural Variant Caller Module
 *
 * Runs SvABA local de novo assembly on BAM alignments for exome target regions.
 * Standardizes output with BGZIP_SORT_INDEX_SV_TUPLE_VCF to produce (*.svaba.sv.vcf.gz).
 */
process RUN_SVABA_ASSEMBLY {
    tag "$sample_id"
    label 'svaba'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fasta_fai
    path target_bed

    output:
    tuple val(sample_id), path("${sample_id}.svaba.raw.sv.vcf"), emit: raw_sv_vcf
    tuple val(sample_id), path("${sample_id}.svaba.indel.vcf"),  emit: indel_vcf, optional: true

    script:
    def is_cram = bam.name.endsWith('.cram')
    """
    #!/bin/bash
    set -euo pipefail

    # 1. Ensure BWA indices exist (SvABA uses internal BWA-MEM during assembly)
    if [ ! -f "${fasta}.bwt" ]; then
        if [ -f "\$(dirname \$(readlink -f ${fasta}))/\$(basename ${fasta}).bwt" ]; then
            REF_DIR=\$(dirname \$(readlink -f ${fasta}))
            REF_BASE=\$(basename ${fasta})
            ln -sf \${REF_DIR}/\${REF_BASE}.* . 2>/dev/null || true
        fi
        if [ ! -f "${fasta}.bwt" ]; then
            echo "[INFO] Indexing ${fasta} with BWA..."
            bwa index ${fasta}
        fi
    fi

    # 2. Alignment Input Handling (BAM vs CRAM)
    ALIGN_IN="${bam}"
    if [ "${is_cram}" = "true" ]; then
        echo "[INFO] Preparing CRAM target regions for SvABA..."
        samtools view -@ ${task.cpus} -b -T ${fasta} -L ${target_bed} ${bam} -o ${sample_id}.target.bam
        samtools index -b ${sample_id}.target.bam
        ALIGN_IN="${sample_id}.target.bam"
    fi

    # 3. Execute SvABA de novo string graph assembly (target-restricted exome mode)
    svaba run \\
        -t \${ALIGN_IN} \\
        -G ${fasta} \\
        -k ${target_bed} \\
        --germline \\
        -a ${sample_id} \\
        -p ${task.cpus}

    if [ "${is_cram}" = "true" ]; then
        rm -f ${sample_id}.target.bam ${sample_id}.target.bam.bai
    fi

    # 4. Extract raw uncompressed SV VCF for standardization via bgzip-sort-index.nf
    if [ -f ${sample_id}.sv.vcf.gz ]; then
        gzip -dc ${sample_id}.sv.vcf.gz > ${sample_id}.svaba.raw.sv.vcf
    elif [ -f ${sample_id}.svaba.sv.vcf.gz ]; then
        gzip -dc ${sample_id}.svaba.sv.vcf.gz > ${sample_id}.svaba.raw.sv.vcf
    elif [ -f ${sample_id}.svaba.sv.vcf ]; then
        mv ${sample_id}.svaba.sv.vcf ${sample_id}.svaba.raw.sv.vcf
    elif [ -f ${sample_id}.sv.vcf ]; then
        mv ${sample_id}.sv.vcf ${sample_id}.svaba.raw.sv.vcf
    else
        touch ${sample_id}.svaba.raw.sv.vcf
    fi

    if [ -f ${sample_id}.indel.vcf.gz ]; then
        gzip -dc ${sample_id}.indel.vcf.gz > ${sample_id}.svaba.indel.vcf
    elif [ -f ${sample_id}.svaba.indel.vcf.gz ]; then
        gzip -dc ${sample_id}.svaba.indel.vcf.gz > ${sample_id}.svaba.indel.vcf
    fi
    """
}

workflow SVABA {
    take:
    bams_ch     // tuple val(sample_id), path(bam), path(bai)
    fasta
    fasta_fai
    target_bed

    main:
    RUN_SVABA_ASSEMBLY(bams_ch, fasta, fasta_fai, target_bed)
    BGZIP_SORT_INDEX_SV_TUPLE_VCF(RUN_SVABA_ASSEMBLY.out.raw_sv_vcf, 'SvABA', 'out_SVABA')

    emit:
    sv_vcf           = BGZIP_SORT_INDEX_SV_TUPLE_VCF.out.sv_vcf
    sorted_vcf       = BGZIP_SORT_INDEX_SV_TUPLE_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_SV_TUPLE_VCF.out.sorted_vcf_index
}

workflow RUN_SVABA {
    take:
    bams_ch
    fasta
    fasta_fai
    target_bed

    main:
    SVABA(bams_ch, fasta, fasta_fai, target_bed)

    emit:
    sv_vcf           = SVABA.out.sv_vcf
    sorted_vcf       = SVABA.out.sorted_vcf
    sorted_vcf_index = SVABA.out.sorted_vcf_index
}
