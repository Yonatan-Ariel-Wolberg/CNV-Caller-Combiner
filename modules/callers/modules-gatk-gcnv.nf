#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_TUPLE_VCF as BGZIP_SORT_INDEX_VCF } from '../common/bgzip-sort-index.nf'
outdir = file(params.outdir, type: 'dir')

// Process to generate ploidy priors
process GENERATE_PLOIDY_PRIORS {
    tag "Ploidy_Priors"
    publishDir "${outdir}/out_GCNV/priors", mode: 'copy', overwrite: true

    input: 
    path fai
    
    output: 
    path "ploidy_priors.tsv", emit: priors
    
    script:
    """
    # Create header with 4 ploidy states (0, 1, 2, 3)
    echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" > ploidy_priors.tsv

    # Default for all contigs: high priority for Ploidy 2 (Autosomes)
    awk -v OFS="\t" '{print \$1, "0.01", "0.01", "0.97", "0.01"}' $fai >> ploidy_priors.tsv

    # Adjust Sex Chromosomes (X and Y) for mixed cohort
    sed -i 's/^chrX\t0.01\t0.01\t0.97\t0.01/chrX\t0.01\t0.49\t0.49\t0.01/' ploidy_priors.tsv
    sed -i 's/^chrY\t0.01\t0.01\t0.97\t0.01/chrY\t0.50\t0.50\t0.00\t0.00/' ploidy_priors.tsv
    """
}

// Process to preprocess intervals
process PREPROCESS_INTERVALS {
    tag "Preprocess_Intervals"
    label 'gatk'
    
    input: 
    path fasta
    path fai
    path dict
    path targets
    
    output: 
    path "preprocessed.interval_list", emit: interval_list
    
    script:
    """
    gatk --java-options "-Xmx4g" PreprocessIntervals \\
        -R $fasta \\
        -L $targets \\
        --bin-length ${params.bin_length} \\
        --padding ${params.padding} \\
        -imr OVERLAPPING_ONLY \\
        -O preprocessed.interval_list
    """
}

// Process to annotate intervals
process ANNOTATE_INTERVALS {
    tag "Annotate_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    path fasta
    path fai
    path dict
    
    output: 
    path "annotated.tsv", emit: annotated_intervals
    
    script:
    """
    gatk --java-options "-Xmx4g" AnnotateIntervals \\
        -L $intervals -R $fasta -imr OVERLAPPING_ONLY -O annotated.tsv
    """
}

// Process to collect read counts
process COLLECT_READ_COUNTS {
    tag "$sample_id"
    label 'gatk'
    
    input: 
    tuple val(sample_id), path(reads), path(index)
    path intervals
    path fasta
    path fai
    path dict

    output: 
    tuple val(sample_id), path("${sample_id}.counts.hdf5"), emit: counts

    script:
    """
    # Resolve symlink to get the real source path of the BAM index
    real_idx=\$(readlink -f $index 2>/dev/null || echo "$index")
    cp -f "\$real_idx" ${reads}.bai 2>/dev/null || true
    touch ${reads}.bai

    gatk --java-options "-Xmx4g -Djava.io.tmpdir=\$PWD" \\
        CollectReadCounts \\
        -L $intervals -R $fasta -I $reads \\
        --format HDF5 -imr OVERLAPPING_ONLY \\
        -O ${sample_id}.counts.hdf5
    """
}

// Process to filter intervals
process FILTER_INTERVALS {
    tag "Filter_Intervals"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path annotated
    path counts
    
    output: 
    path "cohort.filtered.interval_list", emit: filtered_intervals
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" FilterIntervals \\
        -L $intervals \\
        --annotated-intervals $annotated \\
        -I counts_files.list \\
        -imr OVERLAPPING_ONLY \\
        -O cohort.filtered.interval_list
    """
}

// Process to determine ploidy for the cohort
process DETERMINE_PLOIDY_COHORT {
    tag "Determine_Ploidy"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path counts
    path priors
    
    output: 
    path "ploidy-calls", emit: ploidy_calls
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" DetermineGermlineContigPloidy \
        -L $intervals \
        -I counts_files.list \
        --contig-ploidy-priors $priors \
        --output . \
        --output-prefix ploidy \
        -imr OVERLAPPING_ONLY
    """
}

// Process to scatter intervals
process SCATTER_INTERVALS {
    tag "Scatter_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    
    output: 
    path "temp_*/scattered.interval_list", emit: shards
    
    script:
    """
    gatk --java-options "-Xmx4g" IntervalListTools \\
        --INPUT $intervals \\
        --SUBDIVISION_MODE INTERVAL_COUNT \\
        --SCATTER_CONTENT ${params.scatter_count} \\
        --OUTPUT .
    """
}

// Process germline CNV calling for the cohort
process GERMLINE_CNV_CALLER_COHORT {
    tag "${shard_prefix}"
    label 'gatk'
    label 'gpu_or_high_cpu'
    
    input: 
    tuple val(shard_prefix), path(interval_shard)
    path annotated
    path counts
    path ploidy_calls
    
    output: 
    path "${shard_prefix}-calls", emit: calls
    path "${shard_prefix}-model", emit: model
    
    script:
    // Output prefix is derived from interval_shard.parent.name to guarantee unique shard directory names
    def sensitivity_params = "--class-coherence-length 1000.0 --cnv-coherence-length 1000.0 --enable-bias-factors false --interval-psi-scale 1.0E-6 --log-mean-bias-standard-deviation 0.01 --sample-psi-scale 1.0E-6"
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" GermlineCNVCaller \\
        --run-mode COHORT \\
        -L $interval_shard \\
        -I counts_files.list \\
        --contig-ploidy-calls $ploidy_calls \\
        --annotated-intervals $annotated \\
        --output . \\
        --output-prefix ${shard_prefix} \\
        -imr OVERLAPPING_ONLY \\
        $sensitivity_params
    """
}

// Process to post-process CNV calls for individual samples
process POSTPROCESS_CALLS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${outdir}/out_GCNV/vcfs", mode: 'copy', overwrite: true
    
    input: 
    tuple val(index), val(sample_id)
    path model_shards
    path call_shards
    path ploidy_calls
    path dict
    
    output: 
    tuple val(sample_id), path("${sample_id}.genotyped_segments.vcf.gz"), emit: final_vcf
    path "${sample_id}.genotyped_intervals.vcf.gz", emit: intervals_vcf
    path "${sample_id}.denoised_copy_ratios.tsv", emit: denoised_copy_ratios
    
    script:
    // Sort shards by directory name so PostprocessGermlineCNVCalls receives them
    // in the correct genomic order (matching the order used by GermlineCNVCaller).
    def sort_shards = { shards -> (shards instanceof List ? shards : [shards]).sort { a, b -> a.name <=> b.name } }
    def sorted_call_shards  = sort_shards(call_shards)
    def sorted_model_shards = sort_shards(model_shards)
    def call_args  = sorted_call_shards .collect { c -> "--calls-shard-path ${c}" }.join(" \\\n        ")
    def model_args = sorted_model_shards.collect { m -> "--model-shard-path ${m}" }.join(" \\\n        ")
    """
    # Dynamically determine the names of allosomal contigs (X and Y)
    if grep -q "SN:chrX" "${dict}"; then
        allosomal_args="--allosomal-contig chrX --allosomal-contig chrY"
    elif grep -q "SN:X" "${dict}"; then
        allosomal_args="--allosomal-contig X --allosomal-contig Y"
    else
        allosomal_args="--allosomal-contig chrX --allosomal-contig chrY"
    fi

    gatk --java-options "-Xmx8g" PostprocessGermlineCNVCalls \\
        ${call_args} \\
        ${model_args} \\
        --contig-ploidy-calls $ploidy_calls \\
        --sample-index $index \\
        \${allosomal_args} \\
        --sequence-dictionary $dict \\
        --output-denoised-copy-ratios ${sample_id}.denoised_copy_ratios.tsv \\
        --output-genotyped-intervals ${sample_id}.genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${sample_id}.genotyped_segments.vcf.gz
    """
}


// Process to enforce sample-specific sex priors by modifying DetermineGermlineContigPloidy's outputs
process ENFORCE_SEX_METADATA {
    tag "Enforce_Sex_Metadata"
    label 'gatk'
    publishDir "${outdir}/out_GCNV/ploidy-calls-corrected", mode: 'copy', overwrite: true

    input:
    path ploidy_calls
    path sexinfo
    val sorted_samples // list of sample IDs in alphabetical order

    output:
    path "corrected-calls", emit: ploidy_calls

    script:
    def sample_list_str = sorted_samples.join(" ")
    """
    echo "Copying original ploidy calls..."
    cp -r ${ploidy_calls} corrected-calls

    # Check if sexinfo file is set and exists
    if [ ! -f "${sexinfo}" ] || [ "${sexinfo.name}" = "NO_SEXINFO_FILE" ]; then
        echo "No sexinfo file provided or found. Bypassing metadata sex priors enforcement."
        exit 0
    fi

    echo "Using sexinfo file: ${sexinfo}"

    # Build a bash array of sample IDs from the alphabetical list
    samples_array=( ${sample_list_str} )

    # Read sexinfo file and update contig_ploidy.tsv for matching samples
    while read -r sample_id sex || [ -n "\$sample_id" ]; do
        # Trim carriage returns and whitespace
        sample_id=\$(echo "\$sample_id" | tr -d '\\r' | xargs)
        sex=\$(echo "\$sex" | tr -d '\\r' | xargs)

        if [ -z "\$sample_id" ]; then
            continue
        fi

        # Find the 0-based index of the sample_id in the alphabetical array
        sample_idx=-1
        for i in "\${!samples_array[@]}"; do
            if [ "\${samples_array[\$i]}" = "\$sample_id" ]; then
                sample_idx=\$i
                break
            fi
        done

        if [ \$sample_idx -eq -1 ]; then
            echo "Warning: sample \$sample_id in sexinfo not found in cohort."
            continue
        fi

        tsv_file="corrected-calls/sample_\${sample_idx}/contig_ploidy.tsv"
        if [ -f "\$tsv_file" ]; then
            echo "Correcting ploidy calls for sample \$sample_id (index \$sample_idx) to sex \$sex"
            if [ "\$sex" = "M" ]; then
                # Male: X -> 1, Y -> 1
                sed -i -E 's/^(chr)?X\\t.*/\\1X\\t1\\t999/' \$tsv_file
                sed -i -E 's/^(chr)?Y\\t.*/\\1Y\\t1\\t999/' \$tsv_file
            elif [ "\$sex" = "F" ]; then
                # Female: X -> 2, Y -> 0
                sed -i -E 's/^(chr)?X\\t.*/\\1X\\t2\\t999/' \$tsv_file
                sed -i -E 's/^(chr)?Y\\t.*/\\1Y\\t0\\t999/' \$tsv_file
            else
                echo "Warning: unknown sex '\$sex' for sample \$sample_id."
            fi
        else
            echo "Warning: contig_ploidy.tsv not found for sample \$sample_id at \$tsv_file"
        fi
    done < "${sexinfo}"
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow GATK_GCNV {
    take:
    bam_ch      // channel: tuple(sample_id, reads, index)
    fasta_ch    // path: reference FASTA
    fai_ch      // path: FASTA index (.fai)
    dict_ch     // path: sequence dictionary (.dict)
    targets_ch  // path: target intervals BED/interval_list

    main:
    // Step 1: Generate ploidy priors from the reference index
    GENERATE_PLOIDY_PRIORS(fai_ch)

    // Step 2: Preprocess intervals with binning and padding
    PREPROCESS_INTERVALS(fasta_ch, fai_ch, dict_ch, targets_ch)

    // Step 3: Annotate intervals with GC content and mappability
    ANNOTATE_INTERVALS(
        PREPROCESS_INTERVALS.out.interval_list,
        fasta_ch, fai_ch, dict_ch
    )

    // Step 4: Collect read counts per sample
    COLLECT_READ_COUNTS(
        bam_ch,
        PREPROCESS_INTERVALS.out.interval_list,
        fasta_ch, fai_ch, dict_ch
    )

    // Sort HDF5 count files alphabetically by filename so that every downstream
    // step (FilterIntervals, DetermineGermlineContigPloidy, GermlineCNVCaller)
    // receives them in the same deterministic order.  The same order is used
    // when assigning sample indices for PostprocessGermlineCNVCalls.
    sorted_counts_ch = COLLECT_READ_COUNTS.out.counts
        .map { sample_id, hdf5 -> hdf5 }
        .toSortedList { a, b -> a.name <=> b.name }

    // Step 5: Filter intervals based on annotation and read count statistics
    FILTER_INTERVALS(
        PREPROCESS_INTERVALS.out.interval_list,
        ANNOTATE_INTERVALS.out.annotated_intervals,
        sorted_counts_ch
    )

    // Step 6: Determine ploidy for the cohort
    DETERMINE_PLOIDY_COHORT(
        FILTER_INTERVALS.out.filtered_intervals,
        sorted_counts_ch,
        GENERATE_PLOIDY_PRIORS.out.priors
    )

    // Step 6b: Enforce metadata-specified sex priors (Method 1)
    // Check if sexinfo parameter is defined
    def sexinfo_file = params.get('sexinfo', '') ? file(params.sexinfo) : file("NO_SEXINFO_FILE")

    sorted_samples_ch = bam_ch
        .map { sample_id, reads, index -> sample_id }
        .toSortedList()

    ENFORCE_SEX_METADATA(
        DETERMINE_PLOIDY_COHORT.out.ploidy_calls,
        sexinfo_file,
        sorted_samples_ch
    )

    // Step 7: Scatter filtered intervals into shards for parallel processing
    SCATTER_INTERVALS(FILTER_INTERVALS.out.filtered_intervals)

    shards_ch = SCATTER_INTERVALS.out.shards.flatten()
        .map { file -> tuple(file.parent.name, file) }

    // Step 8: Run germline CNV caller in cohort mode for each interval shard
    GERMLINE_CNV_CALLER_COHORT(
        shards_ch,
        ANNOTATE_INTERVALS.out.annotated_intervals,
        sorted_counts_ch,
        ENFORCE_SEX_METADATA.out.ploidy_calls
    )

    // Step 9: Post-process calls into per-sample genotyped VCFs.
    // Sample indices must be assigned in the same alphabetical order used when
    // building counts_files.list for DetermineGermlineContigPloidy.
    indexed_samples_ch = bam_ch
        .map { sample_id, reads, index -> sample_id }
        .toSortedList()
        .flatMap { samples ->
            samples.withIndex().collect { sample, idx -> tuple(idx, sample) }
        }

    POSTPROCESS_CALLS(
        indexed_samples_ch,
        GERMLINE_CNV_CALLER_COHORT.out.model.collect(),
        GERMLINE_CNV_CALLER_COHORT.out.calls.collect(),
        ENFORCE_SEX_METADATA.out.ploidy_calls,
        dict_ch
    )

    // Step 10: Compress, sort, index, and annotate each segments VCF with TOOL=GATK-gCNV.
    // POSTPROCESS_CALLS now emits tuple(sample_id, vcf) so BGZIP_SORT_INDEX_VCF can
    // embed the sample ID in the output filename using the _GCNV convention expected
    // by gather_vcfs() in main.nf.
    BGZIP_SORT_INDEX_VCF(POSTPROCESS_CALLS.out.final_vcf, 'GATK-gCNV', 'out_GCNV')

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    intervals_vcf    = POSTPROCESS_CALLS.out.intervals_vcf
}
