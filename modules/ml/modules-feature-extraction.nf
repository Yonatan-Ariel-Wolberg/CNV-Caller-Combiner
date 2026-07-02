#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESS: EXTRACT FEATURES
// =====================================================================================
//
// Runs bin/feature_extraction.py on a SURVIVOR- or Truvari-merged VCF and the
// per-caller normalised VCFs produced by normalise_cnv_caller_quality_scores.py.
//
// All per-caller VCF inputs and auxiliary annotation files are OPTIONAL.
// Features that require absent files will be NaN in the output TSV rather
// than causing the process to fail.
//
// Input tuple layout (all optional fields default to [] / '' when absent):
//   val  sample_id
//   path merged_vcf        -- SURVIVOR or Truvari merged VCF
//   path collapsed_vcf     -- Truvari collapsed VCF (truvari collapse -c output)
//                             Pass [] when not available (SURVIVOR mode or absent).
//                             Used with MatchId to trace which callers contributed
//                             to each representative merged call in Truvari mode.
//   val  tool_vcfs_str     -- comma-separated caller=vcf_path pairs
//                             (empty string '' if no per-caller VCFs available)
//   val  merger_mode       -- 'survivor' | 'truvari'
//   path bam_file          -- BAM/CRAM for RD-ratio and L2R stats  (or [])
//   path reference_fasta   -- indexed FASTA for GC content          (or [])
//   path bed_file          -- capture target BED for probe counts   (or [])
//   path mappability_file  -- 4-column mappability BED              (or [])
//   path indelible_counts  -- INDELIBLE count TSV                   (or [])

process EXTRACT_FEATURES {
    tag "${sample_id}"
    label 'feature_extraction'
    publishDir "${outdir}/out_FEATURES/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(merged_vcf),
          path(collapsed_vcf),
          val(tool_vcfs_str),
          val(merger_mode),
          path(bam_file),
          path(reference_fasta),
          path(bed_file),
          path(mappability_file),
          path(indelible_counts)

    output:
    tuple val(sample_id), path("${sample_id}_features.tsv"), emit: features_tsv

    script:
    def bam_arg           = bam_file          ? "--bam_file '${bam_file}'"                 : ''
    def fasta_arg         = reference_fasta   ? "--reference_fasta '${reference_fasta}'"   : ''
    def bed_arg           = bed_file          ? "--bed_file '${bed_file}'"                 : ''
    def map_arg           = mappability_file  ? "--mappability_file '${mappability_file}'" : ''
    def indelible_arg     = indelible_counts  ? "--indelible_counts '${indelible_counts}'" : ''
    def tool_vcfs_arg     = tool_vcfs_str     ? "--tool_vcfs '${tool_vcfs_str}'"           : ''
    def collapsed_arg     = collapsed_vcf     ? "--collapsed_vcf '${collapsed_vcf}'"       : ''
    """
    python ${projectDir}/bin/feature_extraction.py \\
        --merged_vcf    '${merged_vcf}'            \\
        --output        '${sample_id}_features.tsv' \\
        --merger_mode   '${merger_mode}'           \\
        --sample_id     '${sample_id}'             \\
        ${collapsed_arg} \\
        ${tool_vcfs_arg} \\
        ${bam_arg}       \\
        ${fasta_arg}     \\
        ${bed_arg}       \\
        ${map_arg}       \\
        ${indelible_arg}
    """
}

// =====================================================================================
// SUB-WORKFLOW
// =====================================================================================

workflow FEATURE_EXTRACTION {
    take:
        // Channel of tuples:
        // [ sample_id, merged_vcf, collapsed_vcf, tool_vcfs_str, merger_mode,
        //   bam_file, reference_fasta, bed_file, mappability_file, indelible_counts ]
        feature_inputs

    main:
        EXTRACT_FEATURES(feature_inputs)

    emit:
        features_tsv = EXTRACT_FEATURES.out.features_tsv
}
