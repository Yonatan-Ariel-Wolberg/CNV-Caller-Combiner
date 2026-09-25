#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// ALLELE FREQUENCY FILTERS MODULE
// Filters merged consensus callsets by cohort allele frequency (AF) thresholds
// to extract rare germline CNV candidate regions (default AF <= 5%).
// Focuses strictly and internally on our dataset/cohort samples.
// It does NOT use external databases (such as ClinVar, gnomAD, 1000G, or DGV).
// External annotations/ACMG classifications are decoupled to modules/annotation/.
// =====================================================================================

process FILTER_RARE_CNVS {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(cohort_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), path('rare_regions.tsv'), emit: rare_data

    script:
    """
    python3 "${projectDir}/bin/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --total-samples ${params.get('total_samples', 90)} \\
        --max-af ${params.get('max_af', 0.05)} \\
        --output-regions "rare_regions.tsv"
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
