#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// COMMON PROCESSES
// Shared processes used by multiple CNV caller modules.
// =====================================================================================

// Annotate a per-sample VCF with TOOL=<tool_annot>, then bgzip-compress,
// coordinate-sort, and tabix-index it.
// Used by callers that do not embed TOOL in their VCF conversion step.
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir { "${params.outdir}/${dir_suffix}/vcfs" }, mode: 'copy', overwrite: true

    input:
    path vcf_file
    val  tool_annot  // annotation string written to INFO/TOOL (e.g. 'CANOES', 'CNVkit')
    val  dir_suffix  // publish directory suffix          (e.g. 'out_CANOES', 'out_CNVKIT')

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.simpleName
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Build a BED annotation file with TOOL=<tool_annot> for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf_file} | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=<tool_annot> in the INFO field
    bcftools annotate \\
        -a ${sample_name}_tool_annot.bed.gz \\
        -c CHROM,FROM,TO,INFO/TOOL \\
        -h extra_header.txt \\
        ${vcf_file} \\
        -O v -o ${sample_name}_annotated.vcf

    # Compress the annotated VCF with bgzip
    bgzip -c ${sample_name}_annotated.vcf > ${sample_name}_annotated.vcf.gz

    # Sort the compressed VCF
    bcftools sort ${sample_name}_annotated.vcf.gz -o ${sorted_gz} -O z

    # Index the sorted VCF
    tabix -p vcf ${sorted_gz}

    # Remove intermediate files
    rm -f ${sample_name}_tool_annot.bed.gz ${sample_name}_tool_annot.bed.gz.tbi \\
          ${sample_name}_annotated.vcf ${sample_name}_annotated.vcf.gz
    """
}

// Normalise CNV quality scores to a common scale.
// Used by all CNV callers except DRAGEN (which has a different publish path).
process NORMALISE_CNV_QUALITY_SCORES {
    tag "${vcf.simpleName}"
    label 'pysam'
    publishDir { "${params.outdir}/${dir_suffix}/vcfs" }, mode: 'copy', overwrite: true

    input:
    path vcf
    val  caller_name  // value passed to --caller (e.g. 'CANOES', 'GATK')
    val  dir_suffix   // publish directory suffix  (e.g. 'out_CANOES', 'out_GCNV')

    output:
    path("*.normalised.vcf.gz"),     emit: normalised_vcf
    path("*.normalised.vcf.gz.tbi"), emit: normalised_vcf_index

    script:
    def sample_name   = vcf.name - '.sorted.vcf.gz'
    def normalised_gz = "${sample_name}.normalised.vcf.gz"
    """
    normalise_cnv_caller_quality_scores.py \\
        --input_vcf ${vcf} \\
        --output_vcf ${sample_name}.normalised.vcf \\
        --caller ${caller_name}
    python3 -c "import pysam; pysam.tabix_index('${sample_name}.normalised.vcf', preset='vcf', force=True)"
    """
}
