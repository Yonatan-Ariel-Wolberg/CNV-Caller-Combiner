#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESSES FOR TRUVARI
// =====================================================================================

process MERGE_VCFS {
    tag "${sample_id}"
    label 'truvari'
    
    input:
    // Takes a sample ID and a list of VCF files (from multiple callers) associated with that sample
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("merged_${sample_id}.vcf.gz"), path("merged_${sample_id}.vcf.gz.tbi"), emit: merged_data

    script:
    """
    # Initialize a string to hold the list of sorted VCFs for concatenation
    sorted_vcfs=""

    for vcf in ${vcfs}; do
        bgzip_vcf="\${vcf}.gz"
        sorted_vcf="\${vcf}.sorted.gz"
        
        # Check if VCF is already bgzipped
        if [ ! -f "\$bgzip_vcf" ]; then
            # 1. Remove records with empty CHROM fields (prevents bcftools errors)
            grep -v "^#" "\$vcf" | awk '\$1 != ""' > cleaned.tmp
            grep "^#" "\$vcf" > header.tmp
            cat header.tmp cleaned.tmp > "cleaned_\$vcf"
        
            # 2. Compress the cleaned VCF
            bgzip -c "cleaned_\$vcf" > "\$bgzip_vcf"
        fi
        
        # Check if the VCF is sorted
        if [ ! -f "\$sorted_vcf" ]; then
            # 3. Sort the compressed VCF
            bcftools sort -o "\$sorted_vcf" -O z "\$bgzip_vcf"
        fi
        
        # 4. Check if the sorted VCF is indexed; if not, index it
        if [ ! -f "\$sorted_vcf.tbi" ]; then
            tabix -p vcf "\$sorted_vcf"
        fi

        # Append to the list of files to merge
        sorted_vcfs="\$sorted_vcfs \$sorted_vcf"
    done

    # 5. Merge all sorted VCF files for this sample into one
    merged_vcf="merged_${sample_id}.vcf.gz"
    bcftools concat -a -O z -o "\$merged_vcf" \$sorted_vcfs
    
    # 6. Index the merged VCF
    tabix -p vcf "\$merged_vcf"
    """
}

process COLLAPSE_VCFS {
    tag "${sample_id}"
    label 'truvari'
    publishDir "${outdir}/out_TRUVARI/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(merged_vcf), path(merged_vcf_tbi)

    output:
    tuple val(sample_id), path("${sample_id}_truvari_merged.vcf"), emit: merged_vcf
    tuple val(sample_id), path("${sample_id}_truvari_collapsed.vcf"), emit: collapsed_vcf

    script:
    """
    # 7. Collapse with Truvari to generate consensus calls
    truvari collapse \\
        --intra \\
        --debug \\
        --pctseq 0 \\
        --pctsize 0.5 \\
        --pctovl 0.5 \\
        --sizemin 50 \\
        --sizemax 2000000 \\
        --fast-cluster \\
        -i "${merged_vcf}" \\
        -o "${sample_id}_truvari_merged.vcf" \\
        -c "${sample_id}_truvari_collapsed.vcf"
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow TRUVARI {
    take:
        grouped_vcfs

    main:
        grouped_vcfs.ifEmpty { error "TRUVARI received no VCF groups. Check that callers have produced VCF output." }
        MERGE_VCFS(grouped_vcfs)
        COLLAPSE_VCFS(MERGE_VCFS.out.merged_data)

    emit:
        merged_vcf = COLLAPSE_VCFS.out.merged_vcf
        collapsed_vcf = COLLAPSE_VCFS.out.collapsed_vcf
}
