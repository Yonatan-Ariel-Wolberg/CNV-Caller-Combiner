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
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Reheader to ensure sample name inside VCF matches sample_name
    echo "${sample_name}" > sample_id.txt
    bcftools reheader -s sample_id.txt ${vcf_file} -o reheadered.vcf

    # Build a BED annotation file with TOOL=<tool_annot> for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' reheadered.vcf | \\
        awk -v tool="${tool_annot}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, tool}' | \\
        sort -k1,1 -k2,2n | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=<tool_annot> in the INFO field
    bcftools annotate \\
        -a ${sample_name}_tool_annot.bed.gz \\
        -c CHROM,FROM,TO,INFO/TOOL \\
        -h extra_header.txt \\
        reheadered.vcf \\
        -O v -o ${sample_name}_annotated.vcf

    # Compress the annotated VCF with bgzip
    bgzip -c ${sample_name}_annotated.vcf > ${sample_name}_annotated.vcf.gz

    # Sort the compressed VCF
    bcftools sort ${sample_name}_annotated.vcf.gz -o ${sorted_gz} -O z

    # Index the sorted VCF
    tabix -p vcf ${sorted_gz}

    # Remove intermediate files
    rm -f ${sample_name}_tool_annot.bed.gz ${sample_name}_tool_annot.bed.gz.tbi \\
          ${sample_name}_annotated.vcf ${sample_name}_annotated.vcf.gz reheadered.vcf sample_id.txt
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
    path("*.normalised*.cnv.vcf.gz"),     emit: normalised_vcf
    path("*.normalised*.cnv.vcf.gz.tbi"), emit: normalised_vcf_index

    script:
    def caller = caller_name.toLowerCase()
    if (caller == 'gatk') caller = 'gcnv'
    def base_id = vcf.name
        .replaceAll(/\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth).*/, '')
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*/, '')
        .replaceAll(/(\.sorted|\.qs-filtered|\.filtered|\.cnv|\.annotated)*\.(vcf|vcf\.gz)$/, '')
        .replaceAll(/\.md$/, '')
    def is_filtered = vcf.name.contains('.qs-filtered') || vcf.name.contains('_filtered') || vcf.name.contains('.filtered')
    def normalised_base = is_filtered
        ? "${base_id}.${caller}.normalised.qs-filtered.cnv"
        : "${base_id}.${caller}.normalised.cnv"
    def normalised_gz   = "${normalised_base}.vcf.gz"
    """
    normalise_cnv_caller_quality_scores.py \\
        --input_vcf ${vcf} \\
        --output_vcf ${normalised_base}.vcf \\
        --caller ${caller_name}
    python3 -c "import pysam; pysam.tabix_index('${normalised_base}.vcf', preset='vcf', force=True)"
    rm -f ${normalised_base}.vcf
    """
}

// =====================================================================================
// PROCESS: COMPUTE_TARGET_GC_ANNOTATION
// Single-pass GATK AnnotateIntervals for GC content shared across callers (CANOES, XHMM, CLAMMS)
// =====================================================================================
process COMPUTE_TARGET_GC_ANNOTATION {
    tag "compute_gc"
    label 'gatk'
    publishDir "${params.outdir}/reference_gc", mode: 'copy', overwrite: true

    input:
    path target_bed
    path ref_fasta
    path ref_fai

    output:
    path "target_gc_annotated.tsv", emit: gc_table
    path "extreme_gc_targets.txt",   emit: extreme_gc, optional: true

    script:
    def mem = Math.max(2, (task.memory.toGiga() * 0.75).toInteger())
    """
    gatk --java-options "-XX:+UseSerialGC -Xms2g -Xmx${mem}g" \\
        AnnotateIntervals \\
        -R ${ref_fasta} \\
        -L ${target_bed} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O target_gc_annotated.tsv

    # Pre-extract extreme GC targets (< 0.1 or > 0.9) for XHMM
    awk -F'\\t' 'NR>1 && (\$4 < 0.1 || \$4 > 0.9) { print \$1":"\$2"-"\$3 }' target_gc_annotated.tsv > extreme_gc_targets.txt || true
    """
}
