#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// WORKFLOW: BUILD MOSDEPTH CONSENSUS TARGET BED FROM SCRATCH
// =============================================================================
// Runs Mosdepth on each QC-passed sample BAM against the baseline Sanger target
// BED to compute per-interval mean coverage, filters for intervals with >=20x
// coverage across >=90% of samples, merges proximate intervals within 1500 bp,
// and generates the Picard/CLAMMS interval_list via GATK BedToIntervalList.
// =============================================================================

params.samplesheet_bams = ""
params.ref              = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.fai              = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
params.dict             = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.dict"
params.target_bed       = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
params.pass_percent     = 90
params.workflow         = "local"
params.outdir           = "mosdepth_consensus_results"

if (!params.samplesheet_bams) {
    error "ERROR: --samplesheet_bams must be specified."
}

include { RUN_MOSDEPTH; BUILD_CONSENSUS; CONVERT_BED_TO_INTERVAL_LIST } from '../../modules/common/modules-mosdepth-target.nf'

workflow {
    ch_bams = Channel.fromPath(params.samplesheet_bams)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sid = row.SampleID ?: row.sample_id ?: row.sample ?: row.ID
            def bam = file(row.BAM ?: row.bam)
            def bai = file(row.BAI ?: row.bai ?: "${row.BAM ?: row.bam}.bai")
            [ sid, bam, bai ]
        }

    RUN_MOSDEPTH(
        ch_bams,
        file(params.ref),
        file(params.fai),
        file(params.target_bed)
    )

    BUILD_CONSENSUS(
        RUN_MOSDEPTH.out.usable_beds.collect(),
        file(params.target_bed)
    )

    CONVERT_BED_TO_INTERVAL_LIST(
        BUILD_CONSENSUS.out.karyotypic_bed,
        file(params.dict)
    )
}
