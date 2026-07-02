#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================================
//  EVALUATE MODULE
//  Assess CNV caller performance (precision and sensitivity) at the probe level.
// ================================================================================

outdir = file(params.outdir, type: 'dir')

// ---------------------------------------------------------------------------------
// Process 1: Convert a single per-sample VCF to a 5-column BED file
// ---------------------------------------------------------------------------------
process VCF_TO_BED {
    tag "${vcf.simpleName}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path vcf

    output:
    path("*.bed"), emit: bed

    script:
    // Strip common suffixes to obtain a clean sample name
    def sample_name = vcf.name
        .replaceAll(/\.(sorted|annotated|normalised)(\.vcf(\.gz)?)?$/, '')
        .replaceAll(/\.vcf(\.gz)?$/, '')
    """
    vcf_to_bed.py --vcf ${vcf} --output ${sample_name}.bed
    """
}

// ---------------------------------------------------------------------------------
// Process 2: Combine all per-sample BED files into one unified call-set BED
// ---------------------------------------------------------------------------------
process COMBINE_BEDS {
    tag "${caller_name}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path beds
    val  caller_name

    output:
    path("${caller_name}_callset.bed"), emit: combined_bed

    script:
    """
    cat ${beds} > ${caller_name}_callset.bed
    """
}

// ---------------------------------------------------------------------------------
// Process 3: Run evaluate_caller_performance.py and write metrics to a text file
// ---------------------------------------------------------------------------------
process EVALUATE_CALLER {
    tag "${caller_name}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path callset_bed
    path truth_bed
    path probes_bed
    val  caller_name

    output:
    path("${caller_name}_performance.txt"), emit: report

    script:
    """
    evaluate_caller_performance.py \\
        --truth_bed   ${truth_bed}   \\
        --callset_bed ${callset_bed} \\
        --probes_bed  ${probes_bed}  \\
        --output      ${caller_name}_performance.txt
    """
}

// ---------------------------------------------------------------------------------
// Process 4: Evaluate ACMG annotation CSV calls (Sensitivity, Precision, Specificity, F-beta, MCC)
// ---------------------------------------------------------------------------------
process EVALUATE_ACMG_CSV {
    tag "${run_name}-${csv.baseName}"
    publishDir "${params.outdir}/evaluate-results/${run_name}", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(csv), path(probes_bed)
    path map_file
    path truth_bed

    output:
    path "*_performance_report.txt", emit: report

    script:
    def base_name = csv.baseName
    """
    evaluate_separated_helper.py \\
        --input-csv "${csv}" \\
        --map-file "${map_file}" \\
        --truth-bed "${truth_bed}" \\
        --probes-bed "${probes_bed}" \\
        --workflow-type "${run_name}" \\
        --output "${run_name}_${base_name}_performance_report.txt"
    """
}

// ---------------------------------------------------------------------------------
// Process 5: Plot performance comparison metrics dynamic visualization
// ---------------------------------------------------------------------------------
process PLOT_EVALUATION_METRICS {
    publishDir "${params.outdir}/evaluate-results", mode: 'copy', overwrite: true

    input:
    path reports

    output:
    path 'performance_metrics_comparison.png', emit: comparison_plot

    script:
    """
    # Structure files into subdirectories so parser detects Control vs Mosdepth
    mkdir -p Control Mosdepth
    for f in *_performance_report.txt; do
        if [[ "\$f" == *"Control"* ]]; then
            cp "\$f" Control/
        elif [[ "\$f" == *"Mosdepth"* ]]; then
            cp "\$f" Mosdepth/
        else
            cp "\$f" Control/
        fi
    done
    
    plot_evaluation_performance_dynamic.py \\
        --input-dir . \\
        --output-png performance_metrics_comparison.png
    """
}

// ---------------------------------------------------------------------------------
// SUB-WORKFLOW: EVALUATE (Probe-level performance reporting)
// ---------------------------------------------------------------------------------
workflow EVALUATE {
    take:
    vcf_ch      // channel of per-sample VCF paths
    truth_bed   // path: truth set BED
    probes_bed  // path: capture target BED
    caller_name // val:  caller label

    main:
    vcf_ch.ifEmpty { error "EVALUATE received no VCF files." }
    VCF_TO_BED(vcf_ch)
    COMBINE_BEDS(VCF_TO_BED.out.bed.collect(), caller_name)
    EVALUATE_CALLER(
        COMBINE_BEDS.out.combined_bed,
        truth_bed,
        probes_bed,
        caller_name
    )

    emit:
    performance_report = EVALUATE_CALLER.out.report
}

// ---------------------------------------------------------------------------------
// SUB-WORKFLOW: CNV_PERFORMANCE_EVALUATION (ACMG-level statistical metrics & plots)
// ---------------------------------------------------------------------------------
workflow CNV_PERFORMANCE_EVALUATION {
    take:
        csv_ch       // Channel: tuple val(run_name), path(csv), path(probes_bed)
        map_file     // Path: Sanger sample ID mapping file
        truth_bed    // Path: Sanger validated truthset BED file

    main:
        EVALUATE_ACMG_CSV(csv_ch, map_file, truth_bed)
        PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())

    emit:
        performance_reports = EVALUATE_ACMG_CSV.out.report
        comparison_plot     = PLOT_EVALUATION_METRICS.out.comparison_plot
}
