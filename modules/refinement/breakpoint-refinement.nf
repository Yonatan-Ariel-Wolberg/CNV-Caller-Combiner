#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * REFINE_CONSENSUS_BREAKPOINTS
 *
 * Post-consensus breakpoint resolution refinement process.
 * Matches intra-sample depth consensus VCF with sample-matched SV VCF
 * (from DRAGEN-SV or SvABA).
 */
process REFINE_CONSENSUS_BREAKPOINTS {
    tag "$sample_id"
    label 'pysam'
    publishDir "${params.outdir}/consensus_refined", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(consensus_vcf), path(consensus_tbi), path(sv_vcf), path(sv_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.refined.vcf.gz"), path("${sample_id}.refined.vcf.gz.tbi"), emit: refined_vcf
    path("${sample_id}_refinement_summary.txt"), emit: summary

    script:
    """
    python3 ${projectDir}/bin/refine_breakpoints.py \\
        --consensus-vcf ${consensus_vcf} \\
        --sv-vcf ${sv_vcf} \\
        --output-vcf ${sample_id}.refined.vcf.gz \\
        --summary ${sample_id}_refinement_summary.txt \\
        --max-flank ${params.get('refinement_max_flank', 50000)} \\
        --max-delta ${params.get('refinement_max_delta', 25000)} \\
        --focal-threshold ${params.get('refinement_focal_threshold', 10000)} \\
        --min-focal-delta ${params.get('refinement_min_focal_delta', 1000)} \\
        --focal-slack-ratio ${params.get('refinement_focal_slack_ratio', 1.0)} \\
        --min-recip-overlap ${params.get('refinement_min_recip_overlap', 0.20)}
    """
}

workflow RUN_BREAKPOINT_REFINEMENT {
    take:
    consensus_vcf_ch // [sample_id, vcf, tbi]
    sv_vcf_ch        // [sample_id, sv_vcf, sv_tbi]

    main:
    paired_ch = consensus_vcf_ch.join(sv_vcf_ch)
    REFINE_CONSENSUS_BREAKPOINTS(paired_ch)

    emit:
    refined_vcf = REFINE_CONSENSUS_BREAKPOINTS.out.refined_vcf
    summary     = REFINE_CONSENSUS_BREAKPOINTS.out.summary
}
