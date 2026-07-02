#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF; NORMALISE_CNV_QUALITY_SCORES } from '../common/modules-common.nf'

// ================================================================================
//  CANOES MODULE
// ================================================================================

// Define parameters for the workflow
params.ref               = ''
params.probes            = ''
params.canoes_batch_size = 100  // BAMs per bedtools multicov job; reduce if hitting memory limits

ref    = file(params.ref    ?: 'NO_FILE', type: 'file')
probes = file(params.probes ?: 'NO_FILE', type: 'file')
outdir = file(params.outdir, type: 'dir')

// Define the list of chromosomes to process
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

// Process to calculate GC content for each chromosome
process CALC_GC_CANOES {
    tag { "${chr}" }
    label 'gatk'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
    
    input:
    each chr
    
    output:
    tuple val("${chr}"), path("${chr}_gc.txt"), emit: chr_gc_content 

    script:
    mem = task.memory.toGiga() - 4
    
    """
    grep ^"${chr}	" ${probes} > ${chr}_probes_sanger.bed

    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        AnnotateIntervals \
        -R ${ref} \
        -L ${chr}_probes_sanger.bed \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${chr}_gc.tmp

    grep ^"${chr}" ${chr}_gc.tmp | awk '{ print \$1":"\$2"-"\$3"\t"\$4 }' > ${chr}_gc.txt
    """
}

// Process to split the BAM list into fixed-size batches for parallel bedtools runs
process BATCH_BAMS {
    tag { 'batch_bams' }
    label 'canoes_R'

    input:
    path(bam_list)

    output:
    path("bam_batch_*.list"), emit: bam_batches

    script:
    def batch_size = params.canoes_batch_size as int
    """
    sort ${bam_list} | split -l ${batch_size} \\
        --numeric-suffixes=000 --suffix-length=3 --additional-suffix=.list - bam_batch_
    """
}

// Process to generate read counts for one BAM batch per chromosome
process GEN_READ_COUNTS {
    tag { "${chr}_${batch_id}" }
    label 'bedtools'

    input:
    tuple val(batch_id), path(batch_list)
    each chr

    output:
    tuple val("${chr}"), path("${chr}_${batch_id}_canoes_reads.txt"), path("${chr}_${batch_id}_sample_list"), emit: chr_reads_cov

    script:
    """
    ## SUBSET THE PROBES FOR THIS CHROMOSOME
    grep ^"${chr}	" ${probes} > ${chr}_probes_sanger.bed

    ## CREATE SAMPLE LIST FOR THIS BATCH
    while read line
    do
        echo \$(basename \$line) | sed 's/\\.bam\$//'
    done < ${batch_list} > ${chr}_${batch_id}_sample_list

    ## RUN BEDTOOLS FOR THIS BATCH
    mapfile -t bams < ${batch_list}
    bedtools multicov -bams "\${bams[@]}" -bed ${chr}_probes_sanger.bed -q 20 > ${chr}_${batch_id}_canoes_reads.txt
    """
}

// Process to merge per-batch read count matrices into a single per-chromosome matrix
process MERGE_READ_COUNTS {
    tag { "${chr}" }
    label 'bedtools'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true

    input:
    tuple val(chr), path(reads_files), path(sample_lists)

    output:
    tuple val(chr), path("${chr}_canoes_reads.txt"), path("${chr}_sample_list"), emit: chr_reads_cov

    script:
    """
    ## Sort batch files to ensure consistent column ordering across batches
    ls ${chr}_*_canoes_reads.txt | sort > reads_order.txt
    ls ${chr}_*_sample_list      | sort > slist_order.txt

    ## Concatenate sample lists in batch order
    while read f; do
        cat "\$f"
    done < slist_order.txt > ${chr}_sample_list

    ## Extract genomic coordinate columns (1-3) from the first batch
    first_reads=\$(head -1 reads_order.txt)
    cut -f1-3 "\$first_reads" > coords.txt

    ## Extract count columns (col 4+) from every batch, then paste side-by-side
    i=0
    while read f; do
        cut -f4- "\$f" > "counts_\${i}.tmp"
        i=\$((i + 1))
    done < reads_order.txt

    paste coords.txt counts_*.tmp > ${chr}_canoes_reads.txt
    rm -f coords.txt counts_*.tmp reads_order.txt slist_order.txt
    """
}

// Process to run the CANOES algorithm
process RUN_CANOES {
    tag { "${chr}" }
    label 'canoes'
    memory '5 GB'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
    
    input:
    tuple val(chr), path(canoes_reads), path(sample_list), path(gc_content)

    output:
    tuple val(chr), path("${chr}_CNVs_pass.csv"), emit: chr_cnvs_pass
    tuple val(chr), path("${chr}_CNVs_fail.csv"), emit: chr_cnvs_fail
    tuple val(chr), path("${chr}_xcnvs.RData"), emit: chr_cnvs_rdata
    tuple val(chr), path("${chr}_CNV_genotype"), emit: chr_cnvs_geno
    tuple val(chr), path("${chr}_CNV_plots"), emit: chr_cnvs_plots
    
    script:
    """
    sed 's/chr//g; /^X/d; /^Y/d' ${canoes_reads} > ${chr}_canoes_reads_new.txt
    sed 's/chr//g; /^X/d; /^Y/d' ${gc_content} > ${chr}_gc_new.txt 
    
    run_canoes.R ${chr}_gc_new.txt ${chr}_canoes_reads_new.txt ${sample_list} ${task.cpus}
    """
}

// Process to filter CNVs generated by CANOES
process FILTER_CANOES_CNVS {
    tag { "ALL_CNVs" }
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    input:
    path(cnvs)
    
    output:
    path("Sample_CNVs.csv"), emit: all_cnvs
    path("Sample_CNVs_filtered.csv"), emit: filtered_cnvs
    
    script:
    """
    grep --no-filename SAMPLE *_CNVs_pass.csv | uniq > Sample_CNVs.csv
    grep -v --no-filename SAMPLE *_CNVs_pass.csv | sort -k1,1 -k3g,3 >> Sample_CNVs.csv
    head -1 Sample_CNVs.csv > Sample_CNVs_filtered.csv
    awk 'NR>1 && (\$10>=80) && (\$10!="NA") && (\$4>=100) { print }' Sample_CNVs.csv >> Sample_CNVs_filtered.csv
    """
}

// Process to convert the filtered CANOES output to VCF format
process CONVERT_CANOES_TO_VCF {
    tag { "CSV_TO_VCF" }
    label 'pysam'
    publishDir "${outdir}/out_CANOES/vcfs", mode: 'copy', overwrite: true

    input:
    path(filtered_csv)
    path(sample_list)
    path(fai)

    output:
    path("*_CANOES_output.vcf"), emit: vcfs

    script:
    """
    canoes_csv_to_vcf.py \\
        --input_file ${filtered_csv} \\
        --sample_file ${sample_list} \\
        --output_dir . \\
        --fai_file ${fai}
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

// Execute the pipeline
workflow CANOES {
    take:
    bam_list_ch   // path: file containing BAM file paths, one per line
    fai_ch        // path: reference FASTA index (.fai)
    chroms_ch     // channel: chromosome names to process

    main:
    // Step 1: Calculate GC content per chromosome
    CALC_GC_CANOES(chroms_ch)

    // Step 2: Split the BAM list into batches to limit per-job BAM count
    BATCH_BAMS(bam_list_ch)
    bam_batches_ch = BATCH_BAMS.out.bam_batches
        .flatten()
        .map { f -> tuple(f.simpleName, f) }  // (batch_id e.g. "bam_batch_000", batch_file)

    // Step 3: Generate read counts per batch per chromosome
    GEN_READ_COUNTS(bam_batches_ch, chroms_ch)

    // Step 4: Merge per-batch count matrices by chromosome
    merged_reads_ch = GEN_READ_COUNTS.out.chr_reads_cov
        .groupTuple(by: 0)  // group by chr: (chr, [reads_files...], [sample_lists...])

    MERGE_READ_COUNTS(merged_reads_ch)

    // Step 5: Join merged read-counts and GC content by chromosome, then run CANOES
    canoes_input_ch = MERGE_READ_COUNTS.out.chr_reads_cov
        .join(CALC_GC_CANOES.out.chr_gc_content, by: 0)

    RUN_CANOES(canoes_input_ch)

    // Step 6: Collect all per-chromosome pass-CNV files and filter
    all_cnvs_ch = RUN_CANOES.out.chr_cnvs_pass
        .map { chr, cnv_file -> cnv_file }
        .collect()

    FILTER_CANOES_CNVS(all_cnvs_ch)

    // Step 7: Reuse sample list from the first MERGE_READ_COUNTS output
    sample_list_ch = MERGE_READ_COUNTS.out.chr_reads_cov
        .ifEmpty { error "CANOES received no merged chromosome read-count outputs. Check chromosome inputs and BAM batches." }
        .first()
        .map { chr, reads, sample_list -> sample_list }

    CONVERT_CANOES_TO_VCF(
        FILTER_CANOES_CNVS.out.filtered_cnvs,
        sample_list_ch,
        fai_ch
    )

    // Step 8: Compress, sort, index, and annotate each VCF with TOOL=CANOES
    BGZIP_SORT_INDEX_VCF(CONVERT_CANOES_TO_VCF.out.vcfs.flatten(), 'CANOES', 'out_CANOES')

    // Step 9: Normalise quality scores to a common scale
    NORMALISE_CNV_QUALITY_SCORES(BGZIP_SORT_INDEX_VCF.out.sorted_vcf.flatten(), 'CANOES', 'out_CANOES')

    emit:
    sorted_vcf           = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index     = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    normalised_vcf       = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf
    normalised_vcf_index = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf_index
}

