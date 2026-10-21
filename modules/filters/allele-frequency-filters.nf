#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// ALLELE FREQUENCY FILTERS MODULE
// Filters merged consensus callsets by cohort allele frequency (INTERNAL_AF)
// and multi-population maximum allele frequency (AF_POPMAX) across relevant ancestral
// sub-populations (gnomAD v4 AFR, 1000G LWK, gnomAD v4 NFE, and gnomAD v4 SAS).
// This ensures variants common in European (~8%) or Indian/South Asian minorities
// or common in Bantu populations are accurately identified as benign polymorphisms.
// =====================================================================================

params.total_samples        = 90
params.max_af               = 0.05
params.max_pop_af           = 0.01
params.pop_tracks           = null // Flexible string: "AFR=/path/afr.bed,NFE=/path/nfe.bed,..."
params.pop_tracks_manifest  = null // Manifest TSV/CSV: subpopulation_name, track_path
params.pop_vcf              = null // Multi-population VCF (e.g. gnomAD-SV)
params.gnomad_afr           = null
params.thousand_g_lwk       = null
params.gnomad_nfe           = null
params.gnomad_sas           = null
params.gnomad_eas           = null
params.gnomad_amr           = null
params.gnomad_fin           = null
params.gnomad_mid           = null
params.gnomad_asj           = null

process FILTER_RARE_CNVS {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(cohort_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), path('rare_regions.tsv'), emit: rare_data
    tuple val(run_name), path('DDD_Cohort_PopMax_Excluded_CNVs.csv'), emit: excluded_data, optional: true

    script:
    def pop_args = []
    if (params.pop_tracks)          pop_args << "--pop-tracks \"${params.pop_tracks}\""
    if (params.pop_tracks_manifest) pop_args << "--pop-tracks-manifest \"${params.pop_tracks_manifest}\""
    if (params.pop_vcf)             pop_args << "--pop-vcf \"${params.pop_vcf}\""
    if (params.gnomad_afr)          pop_args << "--gnomad-afr \"${params.gnomad_afr}\""
    if (params.thousand_g_lwk)      pop_args << "--thousand-g-lwk \"${params.thousand_g_lwk}\""
    if (params.gnomad_nfe)          pop_args << "--gnomad-nfe \"${params.gnomad_nfe}\""
    if (params.gnomad_sas)          pop_args << "--gnomad-sas \"${params.gnomad_sas}\""
    if (params.gnomad_eas)          pop_args << "--gnomad-eas \"${params.gnomad_eas}\""
    if (params.gnomad_amr)          pop_args << "--gnomad-amr \"${params.gnomad_amr}\""
    if (params.gnomad_fin)          pop_args << "--gnomad-fin \"${params.gnomad_fin}\""
    if (params.gnomad_mid)          pop_args << "--gnomad-mid \"${params.gnomad_mid}\""
    if (params.gnomad_asj)          pop_args << "--gnomad-asj \"${params.gnomad_asj}\""
    def pop_args_str = pop_args.join(" ")
    """
    python3 "${projectDir}/bin/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --output-excluded-csv "DDD_Cohort_PopMax_Excluded_CNVs.csv" \\
        --total-samples ${params.get('total_samples', 90)} \\
        --max-af ${params.get('max_af', 0.05)} \\
        --max-pop-af ${params.get('max_pop_af', 0.01)} \\
        --output-regions "rare_regions.tsv" \\
        ${pop_args_str}
    """
}

process FILTER_VCF_BY_RARE_REGIONS {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf), path(rare_regions)

    output:
    tuple val(run_name), path('frequency-filtered.cohort.cnv.vcf.gz'), emit: rare_vcf
    tuple val(run_name), path('frequency-filtered.cohort.cnv.vcf.gz.tbi'), emit: rare_tbi
    path('cohort_rare.vcf')

    script:
    """
    # Sort and compress regions
    sort -k1,1 -k2,2n ${rare_regions} | bgzip -c > regions.bed.gz
    tabix -p bed regions.bed.gz

    # Ensure input VCF is indexed for region querying
    if [ ! -f "${merged_vcf}.tbi" ] && [ ! -f "${merged_vcf}.csi" ]; then
        tabix -f -p vcf "${merged_vcf}"
    fi

    # Filter cohort VCF by rare regions, outputting harmonized bgzipped and tabix-indexed VCF
    bcftools view -R regions.bed.gz ${merged_vcf} -O z -o frequency-filtered.cohort.cnv.vcf.gz
    tabix -f -p vcf frequency-filtered.cohort.cnv.vcf.gz

    # Uncompressed fallback/alias
    bcftools view frequency-filtered.cohort.cnv.vcf.gz > cohort_rare.vcf
    """
}

// Backward-compatibility wrappers for legacy pipelines
workflow filterRareCnvs {
    take: cohort_ch
    main:
        FILTER_RARE_CNVS(cohort_ch)
    emit:
        rare_data = FILTER_RARE_CNVS.out.rare_data
}

workflow filterVcfByRareRegions {
    take: vcf_and_regions
    main:
        FILTER_VCF_BY_RARE_REGIONS(vcf_and_regions)
    emit:
        rare_vcf = FILTER_VCF_BY_RARE_REGIONS.out.rare_vcf
        rare_tbi = FILTER_VCF_BY_RARE_REGIONS.out.rare_tbi
}
