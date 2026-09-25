nextflow.enable.dsl=2

process RUN_MERGE_CONSENSUS_COHORT {
    label 'survivor'

    input:
    path vcfs

    output:
    path 'cohort_merged.raw.vcf', emit: raw_vcf

    script:
    """
    # Create the list of VCFs to merge
    > vcf_list.txt
    for vcf in ${vcfs}; do
        if [[ "\$vcf" == *.gz ]]; then
            uncompressed="\${vcf%.gz}"
            gunzip -c "\$vcf" > "\$uncompressed"
            echo "\$uncompressed" >> vcf_list.txt
        else
            echo "\$vcf" >> vcf_list.txt
        fi
    done

    # Merge cohort-wide using SURVIVOR with configurable inter-sample distance (default 10,000 bp for refined calls)
    SURVIVOR merge vcf_list.txt ${params.get('survivor_inter_distance', 10000)} 1 1 0 0 30 cohort_merged.raw.vcf
    """
}

process BGZIP_INDEX_COHORT {
    label 'bcftools'
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    path raw_vcf

    output:
    path 'cohort.cnv.vcf.gz', emit: merged_vcf
    path 'cohort.cnv.vcf.gz.tbi', emit: merged_tbi
    path 'cohort_merged.vcf'

    script:
    """
    bcftools sort -O z -o cohort.cnv.vcf.gz ${raw_vcf}
    tabix -f -p vcf cohort.cnv.vcf.gz
    cp ${raw_vcf} cohort_merged.vcf
    """
}

workflow MERGE_CONSENSUS_COHORT {
    take:
        vcfs

    main:
        RUN_MERGE_CONSENSUS_COHORT(vcfs)
        BGZIP_INDEX_COHORT(RUN_MERGE_CONSENSUS_COHORT.out.raw_vcf)

    emit:
        merged_vcf = BGZIP_INDEX_COHORT.out.merged_vcf
        merged_tbi = BGZIP_INDEX_COHORT.out.merged_tbi
}

process SPLIT_COHORT_VCF {
    label 'bcftools'
    tag "${sample_id}"
    publishDir "${params.outdir}/split_rare_vcfs", mode: 'copy'

    input:
    tuple val(sample_id), path(cohort_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_rare.vcf"), emit: rare_vcf

    script:
    """
    bcftools view -s ${sample_id} ${cohort_vcf} -o ${sample_id}_rare.vcf || touch ${sample_id}_rare.vcf
    """
}

process SANITIZE_COHORT_VCF {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf)

    output:
    tuple val(run_name), path('cohort_merged_sanitized.vcf'), emit: sanitized_vcf

    script:
    """
    # Replace non-standard genotype strings like "./." or "." with "0/0"
    # This prevents Downstream bcftools query/awk flattening failures
    if [[ "${merged_vcf}" == *.gz ]]; then
        bcftools view "${merged_vcf}" | sed -E 's/\\t\\.\\/\\.:/\\t0\\/0:/g; s/\\t\\.:/\\t0\\/0:/g' > cohort_merged_sanitized.vcf
    else
        sed -E 's/\\t\\.\\/\\.:/\\t0\\/0:/g; s/\\t\\.:/\\t0\\/0:/g' ${merged_vcf} > cohort_merged_sanitized.vcf
    fi
    """
}

process BCFTOOLS_FLATTEN_COHORT {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf)

    output:
    tuple val(run_name), path('cohort_cnvs_by_sample.tsv'), emit: cohort_tsv

    script:
    """
    bcftools query -f '[%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%TY\n]' ${merged_vcf} \
        | awk 'BEGIN{OFS="\t"} \$7 != "NaN" && \$7 != "0/0" && \$7 != "0" && \$7 != "." && \$7 != "" {print \$1, \$2, \$3, \$4, \$5, \$6}' \
        > cohort_cnvs_by_sample.tsv
    """
}

// Backward-compatibility wrappers for legacy pipelines
workflow sanitizeCohortVcf {
    take: vcf_ch
    main:
        SANITIZE_COHORT_VCF(vcf_ch)
    emit:
        sanitized_vcf = SANITIZE_COHORT_VCF.out.sanitized_vcf
}

workflow bcftoolsFlattenCohort {
    take: vcf_ch
    main:
        BCFTOOLS_FLATTEN_COHORT(vcf_ch)
    emit:
        cohort_tsv = BCFTOOLS_FLATTEN_COHORT.out.cohort_tsv
}

include { BUILD_COHORT_FEATURE_DATABASE } from '../ml/modules-cohort-database.nf'
