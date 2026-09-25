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
    label 'bcftools'
    
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
        if [[ "\$vcf" == *.gz ]]; then
            sorted_vcf="\$(basename \$vcf .gz).sorted.gz"
            bcftools sort -o "\$sorted_vcf" -O z "\$vcf"
        else
            bgzip_vcf="\${vcf}.gz"
            sorted_vcf="\${vcf}.sorted.gz"
            
            # Remove records with empty CHROM fields
            grep -v "^#" "\$vcf" | awk '\$1 != ""' > cleaned.tmp
            grep "^#" "\$vcf" > header.tmp
            cat header.tmp cleaned.tmp > "cleaned_\$vcf"
            
            bgzip -c "cleaned_\$vcf" > "\$bgzip_vcf"
            bcftools sort -o "\$sorted_vcf" -O z "\$bgzip_vcf"
        fi
        
        # Rename sample column to sample_id to avoid "Different sample names" error during bcftools concat
        echo "${sample_id}" > sample_name.txt
        bcftools reheader -s sample_name.txt -o "reheadered_\$sorted_vcf" "\$sorted_vcf"
        mv "reheadered_\$sorted_vcf" "\$sorted_vcf"

        tabix -f -p vcf "\$sorted_vcf"
        sorted_vcfs="\$sorted_vcfs \$sorted_vcf"
    done

    # 5. Merge all sorted VCF files for this sample into one
    merged_vcf="merged_${sample_id}.vcf.gz"
    bcftools concat -a -O z -o "\$merged_vcf" \$sorted_vcfs
    
    # Split multi-allelic records (e.g. from CLAMMS) to avoid Truvari collapse error
    bcftools norm -m -any -O z -o "split_\$merged_vcf" "\$merged_vcf"
    mv "split_\$merged_vcf" "\$merged_vcf"

    # 6. Index the merged VCF
    tabix -f -p vcf "\$merged_vcf"
    """
}

process COLLAPSE_VCFS {
    tag "${sample_id}"
    label 'truvari'
    publishDir "${outdir}/out_TRUVARI/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(merged_vcf), path(merged_vcf_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.intra-merge.cnv.vcf.gz"), emit: merged_vcf
    tuple val(sample_id), path("${sample_id}.intra-merge.cnv.vcf.gz.tbi"), emit: merged_tbi
    tuple val(sample_id), path("${sample_id}_truvari_collapsed.vcf"), emit: collapsed_vcf
    path("${sample_id}_truvari_merged.vcf")

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

    python3 ${projectDir}/bin/annotate_intra_merged_vcf.py \\
        --input-vcf "${sample_id}_truvari_merged.vcf" \\
        --output-vcf "${sample_id}.intra-merge.cnv.vcf.gz" \\
        --sample-id "${sample_id}" \\
        --collapsed-vcf "${sample_id}_truvari_collapsed.vcf" \\
        --merger-mode truvari
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
        merged_tbi = COLLAPSE_VCFS.out.merged_tbi
        collapsed_vcf = COLLAPSE_VCFS.out.collapsed_vcf
}
