#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF } from '../common/bgzip-sort-index.nf'

// ================================================================================
//  EXOMEDEPTH MODULE
//  Integrates the ExomeDepth beta-binomial CNV calling algorithm into the pipeline.
// ================================================================================

params.ref                        = ''
params.probes                     = ''
params.fai                        = ''
params.exomedepth_batch_size      = 50    // BAMs per bedtools multicov batch
params.exomedepth_min_bf          = 10.0  // Minimum Bayes Factor for PASS filter
params.exomedepth_min_correlation = 0.97  // Minimum correlation for reference set selection
params.exomedepth_transition_prob = 1e-4  // HMM transition probability
params.outdir                     = './results'

probes = file(params.probes ?: 'NO_FILE', type: 'file')
fai    = file(params.fai    ?: 'NO_FILE', type: 'file')
outdir = file(params.outdir, type: 'dir')

// Process to prepare standard 4-column BED for ExomeDepth targets (chr, start, end, name)
process PREPARE_EXOMEDEPTH_TARGETS {
    tag { "prep_targets" }
    label 'bedtools'
    publishDir "${outdir}/out_EXOMEDEPTH", mode: 'copy', overwrite: true

    input:
    path probes_bed

    output:
    path("exomedepth_targets.bed"), emit: targets_bed

    script:
    """
    # Ensure BED has 4 columns: chrom, start, end, name (fallback name to chrom:start-end)
    awk 'BEGIN{OFS="\\t"} {
        if (NF >= 4) {
            print \$1, \$2, \$3, \$4
        } else {
            print \$1, \$2, \$3, \$1":"\$2"-"\$3
        }
    }' ${probes_bed} | LC_ALL=C sort -k1,1 -k2,2n > exomedepth_targets.bed
    """
}

// Process to split the BAM list into batches to limit per-job memory & open file handles
process BATCH_BAMS_EXOMEDEPTH {
    tag { 'batch_bams' }
    label 'canoes_R'

    input:
    path bam_list

    output:
    path("bam_batch_*.list"), emit: bam_batches

    script:
    def batch_size = params.exomedepth_batch_size as int
    """
    sort ${bam_list} | split -l ${batch_size} \\
        --numeric-suffixes=000 --suffix-length=3 --additional-suffix=.list - bam_batch_
    """
}

// Process to generate read counts for one BAM batch across all target intervals
process GEN_EXOMEDEPTH_COUNTS {
    tag { "${batch_id}" }
    label 'bedtools'

    input:
    tuple val(batch_id), path(batch_list)
    path targets_bed

    output:
    tuple val(batch_id), path("${batch_id}_exomedepth_reads.txt"), path("${batch_id}_sample_list.txt"), emit: batch_counts

    script:
    """
    ## Generate sample list for this batch
    while read -r line; do
        basename "\$line" | sed -E 's/\\.(bam|cram)\$//' | sed 's/\\.md\$//'
    done < ${batch_list} > ${batch_id}_sample_list.txt

    ## Count reads across target intervals using bedtools multicov
    mapfile -t bams < ${batch_list}
    bedtools multicov -bams "\${bams[@]}" -bed ${targets_bed} -q 20 > ${batch_id}_exomedepth_reads.txt
    """
}

// Process to assemble per-batch counts into a unified cohort count matrix
process MERGE_EXOMEDEPTH_COUNTS {
    tag { "merge_counts" }
    label 'bedtools'
    publishDir "${outdir}/out_EXOMEDEPTH", mode: 'copy', overwrite: true

    input:
    path reads_files
    path sample_lists

    output:
    path("exomedepth_cohort_counts.txt"), emit: count_matrix
    path("cohort_sample_list.txt"),       emit: sample_list

    script:
    """
    ls *_exomedepth_reads.txt | sort > reads_order.txt
    ls *_sample_list.txt      | sort > slist_order.txt

    ## Concatenate sample lists in deterministic order
    while read -r f; do
        cat "\$f"
    done < slist_order.txt > cohort_sample_list.txt

    ## Build count matrix header
    header="chromosome\\tstart\\tend\\tname"
    while read -r s; do
        header="\${header}\\t\${s}"
    done < cohort_sample_list.txt

    echo -e "\$header" > exomedepth_cohort_counts.txt

    ## Extract coordinates (cols 1-4) from first batch
    first_reads=\$(head -1 reads_order.txt)
    cut -f1-4 "\$first_reads" > coords.txt

    ## Extract count columns (col 5+) from each batch and paste horizontally
    i=0
    while read -r f; do
        cut -f5- "\$f" > "counts_\${i}.tmp"
        i=\$((i + 1))
    done < reads_order.txt

    paste coords.txt counts_*.tmp >> exomedepth_cohort_counts.txt
    rm -f coords.txt counts_*.tmp reads_order.txt slist_order.txt
    """
}

// Process to run the ExomeDepth algorithm for a single target sample against the cohort matrix
process RUN_EXOMEDEPTH {
    tag { "${sample_id}" }
    label 'exomedepth'
    memory '8 GB'
    cpus 2
    publishDir "${outdir}/out_EXOMEDEPTH/raw_calls", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(count_matrix)

    output:
    tuple val(sample_id), path("${sample_id}_exomedepth_cnvs.csv"), emit: sample_cnvs

    script:
    def min_corr   = params.exomedepth_min_correlation ?: 0.97
    def trans_prob = params.exomedepth_transition_prob ?: 1e-4
    """
    run_exomedepth.R \\
        ${count_matrix} \\
        ${sample_id} \\
        ${sample_id}_exomedepth_cnvs.csv \\
        ${min_corr} \\
        ${trans_prob}
    """
}

// Process to convert ExomeDepth CSV calls to standard VCF format
process CONVERT_EXOMEDEPTH_TO_VCF {
    tag { "${sample_id}" }
    label 'pysam'
    publishDir "${outdir}/out_EXOMEDEPTH/vcfs", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(csv_file)
    path fai_file

    output:
    path("${sample_id}_EXOMEDEPTH_output.vcf"), emit: vcf

    script:
    def min_bf = params.exomedepth_min_bf ?: 10.0
    """
    exomedepth_csv_to_vcf.py \\
        --input_file ${csv_file} \\
        --sample_id ${sample_id} \\
        --output_dir . \\
        --fai_file ${fai_file} \\
        --min_bf ${min_bf}
    """
}

// ================================================================================
//  EXOMEDEPTH SUB-WORKFLOW
// ================================================================================

workflow EXOMEDEPTH {
    take:
    bam_list_ch   // channel: path to file containing BAM paths (one per line)
    probes_bed_ch // channel: path to target intervals BED file
    fai_ch        // channel: path to reference FASTA index (.fai)

    main:
    // Step 1: Format target BED to 4-column BED
    PREPARE_EXOMEDEPTH_TARGETS(probes_bed_ch)

    // Step 2: Split BAM list into manageable batches
    BATCH_BAMS_EXOMEDEPTH(bam_list_ch)
    bam_batches_ch = BATCH_BAMS_EXOMEDEPTH.out.bam_batches
        .flatten()
        .map { f -> tuple(f.simpleName, f) }

    // Step 3: Compute read coverage per batch across targets
    GEN_EXOMEDEPTH_COUNTS(bam_batches_ch, PREPARE_EXOMEDEPTH_TARGETS.out.targets_bed)

    // Step 4: Merge batch counts into unified cohort count matrix
    reads_collected_ch = GEN_EXOMEDEPTH_COUNTS.out.batch_counts
        .map { id, reads, slist -> reads }
        .collect()
    slist_collected_ch = GEN_EXOMEDEPTH_COUNTS.out.batch_counts
        .map { id, reads, slist -> slist }
        .collect()

    MERGE_EXOMEDEPTH_COUNTS(reads_collected_ch, slist_collected_ch)

    // Step 5: Extract sample IDs and pair each sample with the cohort count matrix
    sample_ids_ch = MERGE_EXOMEDEPTH_COUNTS.out.sample_list
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }

    sample_matrix_ch = sample_ids_ch.combine(MERGE_EXOMEDEPTH_COUNTS.out.count_matrix)

    // Step 6: Run ExomeDepth in parallel across all cohort samples
    RUN_EXOMEDEPTH(sample_matrix_ch)

    // Step 7: Convert CSV calls to per-sample VCF
    CONVERT_EXOMEDEPTH_TO_VCF(
        RUN_EXOMEDEPTH.out.sample_cnvs,
        fai_ch
    )

    // Step 8: Compress, coordinate-sort, and index VCFs with TOOL=EXOMEDEPTH
    BGZIP_SORT_INDEX_VCF(
        CONVERT_EXOMEDEPTH_TO_VCF.out.vcf,
        'EXOMEDEPTH',
        'out_EXOMEDEPTH'
    )

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
