#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESSES FOR SURVIVOR
// =====================================================================================

process RUN_SURVIVOR_MERGE {
    tag "${sample_id}"
    label 'survivor'

    input:
    // Takes a sample ID and a list of VCF files (from multiple callers) associated with that sample
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}_survivor_union.vcf"), emit: raw_union_vcf
    tuple val(sample_id), path("${sample_id}_survivor_intersection.vcf"), emit: raw_intersection_vcf
    tuple val(sample_id), path(vcfs), emit: input_vcfs

    script:
    """
    # Create a list file for SURVIVOR
    list_file="${sample_id}_vcfs.list"
    rm -f "\$list_file"  # Remove any existing list file to avoid conflicts

    # Loop through the VCF files to process
    for vcf in ${vcfs}; do
        fname=\$(basename "\$vcf" .gz)
        if [[ "\$vcf" == *.gz ]]; then
            gunzip -c "\$vcf" > "\${fname}.raw"
        else
            cp "\$vcf" "\${fname}.raw"
        fi
        # Ensure diploid GT for SURVIVOR compatibility (converts haploid '1' to '1/1' to prevent SURVIVOR genotype offset bug)
        awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
            split(\$9, format_keys, ":");
            if (format_keys[1] == "GT") {
                split(\$10, fields, ":");
                gt = fields[1];
                if (gt !~ /[\\/|]/) {
                    if (gt == ".") fields[1] = "./.";
                    else fields[1] = gt "/" gt;
                    s = fields[1];
                    for (i=2; i<=length(fields); i++) s = s ":" fields[i];
                    \$10 = s;
                }
            }
            print;
        }' "\${fname}.raw" > "\${fname}"
        rm -f "\${fname}.raw"
        echo "\${fname}" >> "\$list_file"
    done

    # Args: list max_dist min_support use_type use_strand est_dist min_sv_size output
    # use_strand=0: CNV callers do not produce strand information

    # Union merge: report SVs supported by at least 1 caller
    SURVIVOR merge \\
        \$list_file \\
        1000 \\
        1 \\
        1 \\
        0 \\
        0 \\
        30 \\
        "${sample_id}_survivor_union.vcf"

    # Intersection merge: report SVs supported by at least 2 callers
    SURVIVOR merge \\
        \$list_file \\
        1000 \\
        2 \\
        1 \\
        0 \\
        0 \\
        30 \\
        "${sample_id}_survivor_intersection.vcf"
    """
}

process BGZIP_INDEX_SURVIVOR {
    tag "${sample_id}"
    label 'pysam'
    publishDir "${outdir}/out_SURVIVOR/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(raw_union_vcf), path(raw_intersection_vcf), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}.intra-merge.cnv.vcf.gz"), emit: union_vcf
    tuple val(sample_id), path("${sample_id}.intra-merge.cnv.vcf.gz.tbi"), emit: union_tbi
    tuple val(sample_id), path("${sample_id}.intra-merge-intersection.cnv.vcf.gz"), emit: intersection_vcf
    tuple val(sample_id), path("${sample_id}.intra-merge-intersection.cnv.vcf.gz.tbi"), emit: intersection_tbi
    path("${sample_id}_survivor_union.vcf")
    path("${sample_id}_survivor_intersection.vcf")

    script:
    """
    tool_vcfs_str=\$(ls ${vcfs} | tr '\\n' ',' | sed 's/,\$//')
    python3 ${projectDir}/bin/annotate_intra_merged_vcf.py \\
        --input-vcf "${raw_union_vcf}" \\
        --output-vcf "${sample_id}.intra-merge.cnv.vcf.gz" \\
        --sample-id "${sample_id}" \\
        --tool-vcfs "\$tool_vcfs_str" \\
        --merger-mode survivor

    if [ "${raw_union_vcf}" != "${sample_id}_survivor_union.vcf" ]; then
        cp "${raw_union_vcf}" "${sample_id}_survivor_union.vcf"
    fi

    python3 ${projectDir}/bin/annotate_intra_merged_vcf.py \\
        --input-vcf "${raw_intersection_vcf}" \\
        --output-vcf "${sample_id}.intra-merge-intersection.cnv.vcf.gz" \\
        --sample-id "${sample_id}" \\
        --tool-vcfs "\$tool_vcfs_str" \\
        --merger-mode survivor

    if [ "${raw_intersection_vcf}" != "${sample_id}_survivor_intersection.vcf" ]; then
        cp "${raw_intersection_vcf}" "${sample_id}_survivor_intersection.vcf"
    fi
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow SURVIVOR {
    take:
        grouped_vcfs

    main:
        grouped_vcfs.ifEmpty { error "SURVIVOR received no VCF groups. Check that callers have produced VCF output." }
        RUN_SURVIVOR_MERGE(grouped_vcfs)
        def paired_survivor = RUN_SURVIVOR_MERGE.out.raw_union_vcf
            .join(RUN_SURVIVOR_MERGE.out.raw_intersection_vcf)
            .join(RUN_SURVIVOR_MERGE.out.input_vcfs)
        BGZIP_INDEX_SURVIVOR(paired_survivor)

    emit:
        union_vcf        = BGZIP_INDEX_SURVIVOR.out.union_vcf
        union_tbi        = BGZIP_INDEX_SURVIVOR.out.union_tbi
        intersection_vcf = BGZIP_INDEX_SURVIVOR.out.intersection_vcf
}
