#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF } from '../common/bgzip-sort-index.nf'

// Define output directory
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// HELPER FUNCTION
// =====================================================================================
def is_large_run() {
    // Safely fetch test_size, defaulting to -1 if missing
    def size = params.get('test_size', -1).toInteger()
    return (size > 50 || size == -1)
}

// =====================================================================================
// PROCESSES FOR CNVKIT
// =====================================================================================

process GENERATE_ACCESS {
    tag "Access"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    
    input: 
    path fasta
    
    output: 
    path "access.hg38.bed", emit: access_bed
    
    script: 
    """
    cnvkit.py access $fasta -o access.hg38.bed
    """
}

process AUTOBIN {
    tag "Autobin"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    
    input: 
    path fasta
    path targets
    path access
    path refflat
    path bam 
    
    output: 
    path "*.target.bed", emit: target_bed
    path "*.antitarget.bed", emit: antitarget_bed
    
    script: 
    """
    cnvkit.py autobin -t $targets -g $access --fasta $fasta --annotate $refflat --short-names $bam
    """
}

process COVERAGE {
    tag "$sample_id"
    label 'cnvkit'
    memory { is_large_run() ? '16 GB' : '8 GB' }
    
    input: 
    tuple val(sample_id), path(bam), path(bai)
    path target_bed
    path antitarget_bed
    
    output: 
    tuple val(sample_id), path("${sample_id}.targetcoverage.cnn"), emit: target_cov
    tuple val(sample_id), path("${sample_id}.antitargetcoverage.cnn"), emit: antitarget_cov
    
    script:
    """
    # Resolve the symlink to get the real source path of the index file
    real_bai=\$(readlink -f $bai)
    # Remove the staged symlink
    rm -f $bai
    # Copy the real index file into the local directory
    cp "\$real_bai" $bai
    # Touch the local index file to update its timestamp so it is newer than the BAM file
    touch $bai

    # Run target and antitarget coverage collection in parallel, allocating 2 threads to each
    cnvkit.py coverage $bam $target_bed -o ${sample_id}.targetcoverage.cnn -p 2 &
    cnvkit.py coverage $bam $antitarget_bed -o ${sample_id}.antitargetcoverage.cnn -p 2 &
    wait
    """
}

process CREATE_POOLED_REFERENCE {
    tag "Pooled_Ref"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    memory { is_large_run() ? '64 GB' : '24 GB' }
    
    input: 
    path fasta
    path covs 
    
    output: 
    path "pooled_reference.cnn", emit: ref_cnn
    
    script: 
    """
    cnvkit.py reference *.cnn --fasta $fasta -o pooled_reference.cnn
    """
}

process CALL_CNV {
    tag "$sample_id"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/calls/${sample_id}", mode: 'copy', overwrite: true
    
    input: 
    tuple val(sample_id), path(t_cov), path(t_anticov)
    path reference
    
    output: 
    tuple val(sample_id), path("${sample_id}.cnr"), path("${sample_id}.cns"), path("${sample_id}.call.cns"), emit: results
    path "*.pdf"
    
    script:
    """
    cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
    cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns -p ${task.cpus}
    cnvkit.py call ${sample_id}.cns -o ${sample_id}.call.cns
    cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-scatter.pdf
    cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-diagram.pdf
    """
}

process EXPORT_RESULTS {
    tag "$sample_id"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/vcfs", mode: 'copy', overwrite: true, pattern: "*.bed"
    
    input: 
    tuple val(sample_id), path(cnr), path(cns), path(call_cns)
    path fasta_fai
    
    output: 
    path("${sample_id}_CNVKIT_output.vcf"), emit: vcf
    path("${sample_id}_calls.bed"), emit: bed
    
    script:
    """
    # Export raw CNVkit results
    cnvkit.py export vcf $call_cns -i $sample_id -o temp.vcf
    cnvkit.py export bed $call_cns -i $sample_id -o ${sample_id}_calls.bed

    # Generate contig headers from the FAI file
    awk '{printf("##contig=<ID=%s,length=%s>\\n", \$1, \$2)}' $fasta_fai > contigs.hdr

    # Inject the contig headers right after the ##fileformat line in the VCF
    sed -e '/##fileformat/r contigs.hdr' temp.vcf > ${sample_id}_CNVKIT_output.vcf
    rm -f temp.vcf contigs.hdr
    """
}

process FILTER_CNVKIT_CNVS {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_CNVKIT/vcfs", mode: 'copy', overwrite: true, pattern: "*.vcf"
    
    input:
    path vcf_file
    val q_threshold
    
    output:
    path("${vcf_file.name}"), emit: filtered_vcf
    
    script:
    def sample_name = vcf_file.simpleName - '_CNVKIT_output'
    """
    bcftools filter -i "(FMT/GQ >= ${q_threshold} || FMT/CNQ >= ${q_threshold}) && ((INFO/SVTYPE == 'DEL' && INFO/SVLEN >= -5000000) || (INFO/SVTYPE == 'DUP' && INFO/SVLEN <= 5000000))" ${vcf_file} -o temp.vcf
    mv temp.vcf ${sample_name}_CNVKIT_output.vcf
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow CNVKIT {
    take:
    bam_ch              // channel: tuple(sample_id, bam, bai)
    fasta_ch            // path: reference FASTA
    targets_ch          // path: target BED file
    refflat_ch          // path: refFlat annotation file
    bam_for_autobin_ch  // path: single BAM file used for autobin bin size estimation

    main:
    // Step 1: Generate accessible regions from the reference genome
    GENERATE_ACCESS(fasta_ch)

    // Step 2: Auto-bin targets and antitargets
    AUTOBIN(fasta_ch, targets_ch, GENERATE_ACCESS.out.access_bed, refflat_ch, bam_for_autobin_ch)

    // Step 3: Calculate target and antitarget coverage per sample
    COVERAGE(bam_ch, AUTOBIN.out.target_bed, AUTOBIN.out.antitarget_bed)

    // Step 4: Create pooled reference from all sample coverages
    all_covs_ch = COVERAGE.out.target_cov
        .mix(COVERAGE.out.antitarget_cov)
        .map { sample_id, cov -> cov }
        .collect()

    CREATE_POOLED_REFERENCE(fasta_ch, all_covs_ch)

    // Step 5: Call CNVs per sample using the pooled reference
    cnv_input_ch = COVERAGE.out.target_cov
        .join(COVERAGE.out.antitarget_cov, by: 0)

    CALL_CNV(cnv_input_ch, CREATE_POOLED_REFERENCE.out.ref_cnn)

    // Step 6: Export CNV results as VCF and BED per sample
    fasta_fai_ch = file(fasta_ch.toString() + '.fai')
    EXPORT_RESULTS(CALL_CNV.out.results, fasta_fai_ch)

    // Step 7: Compress, sort, index, and annotate each VCF with TOOL=CNVkit (quality filtering is handled centrally in quality-filters.nf)
    BGZIP_SORT_INDEX_VCF(EXPORT_RESULTS.out.vcf, 'CNVkit', 'out_CNVKIT')

    emit:
    target_cov       = COVERAGE.out.target_cov
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    bed              = EXPORT_RESULTS.out.bed
}
