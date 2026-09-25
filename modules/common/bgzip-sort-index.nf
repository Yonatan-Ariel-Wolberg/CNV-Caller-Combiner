#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE: bgzip-sort-index.nf
// Centralized VCF standardization for Raw CNV and Raw SV Callsets:
//   - Sample ID sanitation and VCF reheadering
//   - INFO/TOOL tag injection (e.g. TOOL=CANOES, TOOL=GATK-gCNV, TOOL=SvABA, TOOL=DRAGEN)
//   - Coordinate sorting (bcftools sort)
//   - Bgzip compression (*.cnv.vcf.gz or *.sv.vcf.gz)
//   - Tabix indexing (*.tbi)
// =====================================================================================

// 1. Process for single CNV VCF file input
// Used by CANOES, CLAMMS, XHMM, CNVkit, ExomeDepth, and DRAGEN CNV
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir { "${params.outdir}/${dir_suffix}/vcfs" }, mode: 'copy', overwrite: true

    input:
    path vcf_file
    val  tool_annot  // e.g. 'CANOES', 'CLAMMS', 'XHMM', 'CNVkit', 'DRAGEN', 'EXOMEDEPTH'
    val  dir_suffix  // e.g. 'out_CANOES', 'out_CNVKIT', 'out_DRAGEN'

    output:
    path("*.cnv.vcf.gz"),     emit: sorted_vcf
    path("*.cnv.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def caller = dir_suffix.replaceAll(/^out_/, '').toLowerCase()
    if (caller == 'gatk') caller = 'gcnv'
    def sample_name = vcf_file.name
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*/, '')
        .replaceAll(/\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth).*/, '')
        .replaceAll(/(\.sorted|\.filtered|\.cnv|\.annotated)*\.(vcf|vcf\.gz)$/, '')
        .replaceAll(/\.vcf$/, '')
        .replaceAll(/\.md$/, '')
    def sorted_gz = "${sample_name}.${caller}.cnv.vcf.gz"
    """
    #!/bin/bash
    set -euo pipefail

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    echo "${sample_name}" > sample_id.txt
    bcftools reheader -s sample_id.txt "${vcf_file}" -o reheadered.vcf

    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' reheadered.vcf | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    bcftools annotate \\
        -a ${sample_name}_tool_annot.bed.gz \\
        -c CHROM,FROM,TO,INFO/TOOL \\
        -h extra_header.txt \\
        reheadered.vcf \\
        -O v -o ${sample_name}_annotated.vcf

    bgzip -c ${sample_name}_annotated.vcf > ${sample_name}_annotated.vcf.gz
    bcftools sort ${sample_name}_annotated.vcf.gz -o ${sorted_gz} -O z
    tabix -p vcf ${sorted_gz}

    rm -f ${sample_name}_tool_annot.bed.gz ${sample_name}_tool_annot.bed.gz.tbi \\
          ${sample_name}_annotated.vcf ${sample_name}_annotated.vcf.gz reheadered.vcf sample_id.txt
    """
}

// 2. Process for tuple CNV VCF input: tuple(sample_id, vcf_file)
// Specifically used by callers like GATK gCNV whose output files lack sample IDs in their filenames
process BGZIP_SORT_INDEX_TUPLE_VCF {
    tag "${sample_id}"
    label 'bcftools'
    publishDir { "${params.outdir}/${dir_suffix}/vcfs" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf_file)
    val tool_annot  // e.g. 'GATK-gCNV'
    val dir_suffix  // e.g. 'out_GCNV'

    output:
    path("*.cnv.vcf.gz"),                                              emit: sorted_vcf
    path("*.cnv.vcf.gz.tbi"),                                          emit: sorted_vcf_index
    tuple val(clean_sample_id), path("*.cnv.vcf.gz"), path("*.cnv.vcf.gz.tbi"), emit: sorted_vcf_tuple

    script:
    def caller = dir_suffix.replaceAll(/^out_/, '').toLowerCase()
    if (caller == 'gatk') caller = 'gcnv'
    clean_sample_id = sample_id.toString().replaceAll(/\.md$/, '')
    def sorted_gz = "${clean_sample_id}.${caller}.cnv.vcf.gz"
    """
    #!/bin/bash
    set -euo pipefail

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    echo "${clean_sample_id}" > sample_id.txt
    bcftools reheader -s sample_id.txt "${vcf_file}" -o reheadered.vcf

    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' reheadered.vcf | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${clean_sample_id}_tool_annot.bed.gz
    tabix -p bed ${clean_sample_id}_tool_annot.bed.gz

    bcftools annotate \\
        -a ${clean_sample_id}_tool_annot.bed.gz \\
        -c CHROM,FROM,TO,INFO/TOOL \\
        -h extra_header.txt \\
        reheadered.vcf \\
        -O v -o ${clean_sample_id}_annotated.vcf

    bgzip -c ${clean_sample_id}_annotated.vcf > ${clean_sample_id}_annotated.vcf.gz
    bcftools sort ${clean_sample_id}_annotated.vcf.gz -o ${sorted_gz} -O z
    tabix -p vcf ${sorted_gz}

    rm -f ${clean_sample_id}_tool_annot.bed.gz ${clean_sample_id}_tool_annot.bed.gz.tbi \\
          ${clean_sample_id}_annotated.vcf ${clean_sample_id}_annotated.vcf.gz reheadered.vcf sample_id.txt
    """
}

// 3. Process for SV tuple input: tuple(sample_id, vcf_file)
// Specifically used by SV callers like SvABA
process BGZIP_SORT_INDEX_SV_TUPLE_VCF {
    tag "${sample_id}"
    label 'bcftools'
    publishDir { "${params.outdir}/${dir_suffix}" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf_file)
    val tool_annot  // e.g. 'SvABA'
    val dir_suffix  // e.g. 'out_SVABA'

    output:
    tuple val(clean_sample_id), path("*.sv.vcf.gz"), path("*.sv.vcf.gz.tbi"), emit: sv_vcf
    path("*.sv.vcf.gz"),                                                      emit: sorted_vcf
    path("*.sv.vcf.gz.tbi"),                                                  emit: sorted_vcf_index

    script:
    def caller = dir_suffix.replaceAll(/^out_/, '').toLowerCase()
    clean_sample_id = sample_id.toString().replaceAll(/\.md$/, '')
    def sorted_gz = "${clean_sample_id}.${caller}.sv.vcf.gz"
    """
    #!/bin/bash
    set -euo pipefail

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    echo "${clean_sample_id}" > sample_id.txt
    bcftools reheader -s sample_id.txt "${vcf_file}" -o reheadered.vcf

    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' reheadered.vcf 2>/dev/null | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${clean_sample_id}_tool_annot.bed.gz || true

    if [ -s ${clean_sample_id}_tool_annot.bed.gz ]; then
        tabix -p bed ${clean_sample_id}_tool_annot.bed.gz
        bcftools annotate \\
            -a ${clean_sample_id}_tool_annot.bed.gz \\
            -c CHROM,FROM,TO,INFO/TOOL \\
            -h extra_header.txt \\
            reheadered.vcf \\
            -O z -o ${sorted_gz}
    else
        bgzip -c reheadered.vcf > ${sorted_gz}
    fi

    bcftools sort ${sorted_gz} -o ${sorted_gz}.tmp -O z
    mv ${sorted_gz}.tmp ${sorted_gz}
    tabix -p vcf ${sorted_gz}

    rm -f ${clean_sample_id}_tool_annot.bed.gz ${clean_sample_id}_tool_annot.bed.gz.tbi \\
          reheadered.vcf sample_id.txt
    """
}

// 4. Process for single SV VCF file input
// Specifically used for standalone or ingested SV VCFs (e.g. DRAGEN SV)
process BGZIP_SORT_INDEX_SV_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir { "${params.outdir}/${dir_suffix}" }, mode: 'copy', overwrite: true

    input:
    path vcf_file
    val  tool_annot  // e.g. 'DRAGEN'
    val  dir_suffix  // e.g. 'out_DRAGEN/sv_vcfs'

    output:
    path("*.sv.vcf.gz"),                                                      emit: sorted_vcf
    path("*.sv.vcf.gz.tbi"),                                                  emit: sorted_vcf_index
    tuple val(sample_name), path("*.sv.vcf.gz"), path("*.sv.vcf.gz.tbi"),     emit: sv_vcf_tuple

    script:
    def caller = dir_suffix.replaceAll(/^out_/, '').replaceAll(/\/.*$/, '').toLowerCase()
    sample_name = vcf_file.name
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH|SVABA).*/, '')
        .replaceAll(/\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth|svaba).*/, '')
        .replaceAll(/(\.sorted|\.filtered|\.cnv|\.sv|\.annotated)*\.(vcf|vcf\.gz)$/, '')
        .replaceAll(/\.vcf$/, '')
        .replaceAll(/\.md$/, '')
    def sorted_gz = "${sample_name}.${caller}.sv.vcf.gz"
    """
    #!/bin/bash
    set -euo pipefail

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    echo "${sample_name}" > sample_id.txt
    bcftools reheader -s sample_id.txt "${vcf_file}" -o reheadered.vcf

    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' reheadered.vcf 2>/dev/null | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz || true

    if [ -s ${sample_name}_tool_annot.bed.gz ]; then
        tabix -p bed ${sample_name}_tool_annot.bed.gz
        bcftools annotate \\
            -a ${sample_name}_tool_annot.bed.gz \\
            -c CHROM,FROM,TO,INFO/TOOL \\
            -h extra_header.txt \\
            reheadered.vcf \\
            -O z -o ${sorted_gz}
    else
        bgzip -c reheadered.vcf > ${sorted_gz}
    fi

    bcftools sort ${sorted_gz} -o ${sorted_gz}.tmp -O z
    mv ${sorted_gz}.tmp ${sorted_gz}
    tabix -p vcf ${sorted_gz}

    rm -f ${sample_name}_tool_annot.bed.gz ${sample_name}_tool_annot.bed.gz.tbi \\
          reheadered.vcf sample_id.txt
    """
}
