#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================================
//  EVALUATE SCORE SEPARATED MODULE
//  Assess CNV caller performance on AnnotSV ACMG Score-Separated CSV files
//  at the probe level against a truth set and capture target probes.
// ================================================================================

// 1. Declare all parameters with default values
params.control_csv_glob = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/control/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-3264cf37-6594-4c0c-b283-6beecac654cf/ddduk-90-run-merged/DDD_Cohort_Rare_CNVs_ACMG*.csv'
params.mosdepth_csv_glob = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/consensus/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-8dc80f47-9ded-4d22-9d9b-938dfb5bcabd/ddduk-90-mosdepth-consensus-run-merged/DDD_Cohort_Rare_CNVs_ACMG*.csv'
params.map_file = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/DDDUK-MAP/EGAF-EGAN-ddduk.map'
params.truth_bed = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed'
params.control_probes_bed = '/dataG/ddd/data/resources/canoes/probes_sanger.bed'
params.mosdepth_probes_bed = '/home/ywolberg/dragen-germline-enrichment-v4-4-6/multiqc-filter-mosdepth/ddduk_karyotypic_consensus.bed'
params.outdir = '/home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/evaluate-results'

// 2. Validate input and exit gracefully if requirements are not met
if (!params.map_file || !file(params.map_file).exists()) {
    log.error "Missing or invalid parameter: params.map_file = ${params.map_file}"
    exit 1
}
if (!params.truth_bed || !file(params.truth_bed).exists()) {
    log.error "Missing or invalid parameter: params.truth_bed = ${params.truth_bed}"
    exit 1
}
if (!params.control_probes_bed || !file(params.control_probes_bed).exists()) {
    log.error "Missing or invalid parameter: params.control_probes_bed = ${params.control_probes_bed}"
    exit 1
}
if (!params.mosdepth_probes_bed || !file(params.mosdepth_probes_bed).exists()) {
    log.error "Missing or invalid parameter: params.mosdepth_probes_bed = ${params.mosdepth_probes_bed}"
    exit 1
}

// ---------------------------------------------------------------------------------
// Process: Evaluate performance on a single ACMG score CSV
// ---------------------------------------------------------------------------------
process EVALUATE_ACMG_CSV {
    tag "${workflow_type}-${csv.baseName}"
    publishDir "${params.outdir}/${workflow_type}", mode: 'copy', overwrite: true

    input:
    tuple path(csv), val(workflow_type), path(probes_bed)
    path map_file
    path truth_bed

    output:
    path "*_performance_report.txt", emit: report

    script:
    def base_name = csv.baseName
    """
    #!/bin/bash
    set -euo pipefail
    python3 "${baseDir}/evaluate_separated_helper.py" \\
        --input-csv "${csv}" \\
        --map-file "${map_file}" \\
        --truth-bed "${truth_bed}" \\
        --probes-bed "${probes_bed}" \\
        --workflow-type "${workflow_type}" \\
        --output "${base_name}_performance_report.txt"
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    touch "${csv.baseName}_performance_report.txt"
    """
}

// ---------------------------------------------------------------------------------
// Main workflow execution
// ---------------------------------------------------------------------------------
workflow {
    // Collect all Control CSVs that match the glob pattern and pair with metadata
    control_csvs_ch = Channel.fromPath(params.control_csv_glob)
        .map { csv -> tuple(csv, 'Control', file(params.control_probes_bed)) }

    // Collect all Mosdepth CSVs that match the glob pattern and pair with metadata
    mosdepth_csvs_ch = Channel.fromPath(params.mosdepth_csv_glob)
        .map { csv -> tuple(csv, 'Mosdepth', file(params.mosdepth_probes_bed)) }

    // Merge both channels to process them in parallel
    all_csvs_ch = control_csvs_ch.concat(mosdepth_csvs_ch)

    // Execute evaluation
    EVALUATE_ACMG_CSV(all_csvs_ch, file(params.map_file), file(params.truth_bed))
}
