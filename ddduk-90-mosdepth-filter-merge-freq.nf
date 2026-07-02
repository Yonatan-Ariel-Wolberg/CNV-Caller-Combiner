#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. Declare all parameters with default values to eliminate undefined warnings
params.input_glob = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/consensus/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-8dc80f47-9ded-4d22-9d9b-938dfb5bcabd/EGAF*/*.cnv.vcf.gz'
params.outdir = '/home/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/consensus/dragen_germline_enrichment_v4-4-6_ddd-uk-batch-8dc80f47-9ded-4d22-9d9b-938dfb5bcabd/ddduk-90-mosdepth-consensus-run-merged'
params.total_samples = 90
params.max_af = 0.05

// 2. Validate input and exit gracefully if requirements are not met
if (!params.input_glob) {
    log.error "Missing required parameter: params.input_glob"
    exit 1
}
if (!params.outdir) {
    log.error "Missing required parameter: params.outdir"
    exit 1
}

// 3. Gracefully check if input files exist before starting
def files = file(params.input_glob)
if (files instanceof List) {
    if (files.isEmpty()) {
        log.error "ERROR: No VCF files found matching pattern: ${params.input_glob}"
        exit 1
    }
} else if (files == null || !files.exists()) {
    log.error "ERROR: No VCF files found matching pattern: ${params.input_glob}"
    exit 1
}

// 4. Expose variables for backward compatibility within the script
def input_glob = params.input_glob
def outdir = params.outdir
def total_samples = params.total_samples as Integer
def max_af = params.max_af as Double

process bcftoolsFilterVcfs {
    label 'bcftools'
    tag 'bcftools-filter'
    publishDir outdir, mode: 'copy', overwrite: true

    input:
    path vcfs

    output:
    path 'vcf_list.txt', emit: vcf_list
    path '*_filtered.cnv.vcf', emit: filtered_vcfs

    script:
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.cnv.vcf.gz; do
        base=\$(basename "\$vcf" .cnv.vcf.gz)
        filtered_vcf="\${base}_filtered.cnv.vcf"
        bcftools view -i 'FILTER="PASS" && QUAL>=15 && FORMAT/BC>=3' \\
            -O v -o "\$filtered_vcf" "\$vcf"
        echo "\$filtered_vcf" >> vcf_list.txt
    done
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    printf "sample\n" > vcf_list.txt
    touch sample_filtered.cnv.vcf
    """
}

process SURVIVORMergeCohort {
    label 'survivor'
    tag 'survivor-merge'
    publishDir outdir, mode: 'copy', overwrite: true

    input:
    path vcf_list
    path filtered_vcfs

    output:
    path 'cohort_merged.vcf', emit: merged_vcf

    script:
    """
    #!/bin/bash
    set -euo pipefail
    SURVIVOR merge vcf_list.txt 1000 1 1 0 0 0 cohort_merged.vcf
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    printf "##fileformat=VCFv4.2\n" > cohort_merged.vcf
    """
}

process bcftoolsFlattenCohort {
    label 'bcftools'
    tag 'bcftools-flatten'
    publishDir outdir, mode: 'copy', overwrite: true

    input:
    path merged_vcf

    output:
    path 'cohort_cnvs_by_sample.tsv', emit: cohort_tsv

    script:
    """
    #!/bin/bash
    set -euo pipefail
    bcftools query -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' cohort_merged.vcf \\
        | awk '\$6 != "0/0" && \$6 != "./." && \$6 != "0"' \\
        > cohort_cnvs_by_sample.tsv
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    printf "sample\tcnv\n" > cohort_cnvs_by_sample.tsv
    """
}

process filterRareCnvs {
    tag 'filter-rare-cnvs'
    publishDir outdir, mode: 'copy', overwrite: false

    input:
    path cohort_tsv

    output:
    path 'DDD_Cohort_Rare_CNVs.csv', emit: rare_cnvs
    path 'rare_regions.tsv', emit: rare_regions

    script:
    """
    #!/bin/bash
    set -euo pipefail
    python3 "${baseDir}/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --total-samples ${total_samples} \\
        --max-af ${max_af} \\
        --output-regions "rare_regions.tsv"
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    printf "sample\n" > DDD_Cohort_Rare_CNVs.csv
    touch rare_regions.tsv
    """
}

process filterVcfByRareRegions {
    label 'bcftools'
    tag 'filter-vcf-regions'
    publishDir outdir, mode: 'copy', overwrite: false

    input:
    path merged_vcf
    path rare_regions

    output:
    path 'cohort_rare.vcf', emit: rare_vcf

    script:
    """
    #!/bin/bash
    set -euo pipefail
    bcftools view -T "${rare_regions}" "${merged_vcf}" -o cohort_rare.vcf
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    touch cohort_rare.vcf
    """
}

process runAnnotSV {
    tag 'annotsv'
    publishDir outdir, mode: 'copy', overwrite: false

    input:
    path rare_vcf

    output:
    path 'cohort_rare.annotated.tsv', emit: annotated_tsv

    script:
    """
    #!/bin/bash
    set -euo pipefail
    # Run AnnotSV on the rare VCF using the host-installed version which has access to the full database
    AnnotSV \\
        -SVinputFile "${rare_vcf}" \\
        -outputDir . \\
        -outputFile cohort_rare.annotated \\
        -genomeBuild GRCh38 \\
        -overwrite 0
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    touch cohort_rare.annotated.tsv
    """
}

process addACMGAnnotations {
    tag 'add-acmg-annotations'
    publishDir outdir, mode: 'copy', overwrite: true

    input:
    path rare_cnvs_csv
    path annotated_tsv

    output:
    path 'DDD_Cohort_Rare_CNVs.csv'
    path 'DDD_Cohort_Rare_CNVs_*.csv', optional: true

    script:
    """
    #!/bin/bash
    set -euo pipefail
    python3 "${baseDir}/add-acmg-annotations.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv"
    """
    stub:
    """
    #!/bin/bash
    set -euo pipefail
    printf "sample\n" > DDD_Cohort_Rare_CNVs.csv
    """
}

workflow {
    vcfs_ch = Channel.fromPath(input_glob, checkIfExists: true).collect()
    bcftoolsFilterVcfs(vcfs_ch)
    SURVIVORMergeCohort(bcftoolsFilterVcfs.out.vcf_list, bcftoolsFilterVcfs.out.filtered_vcfs)
    bcftoolsFlattenCohort(SURVIVORMergeCohort.out.merged_vcf)
    filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
    filterVcfByRareRegions(SURVIVORMergeCohort.out.merged_vcf, filterRareCnvs.out.rare_regions)
    runAnnotSV(filterVcfByRareRegions.out.rare_vcf)
    addACMGAnnotations(filterRareCnvs.out.rare_cnvs, runAnnotSV.out.annotated_tsv)
}
