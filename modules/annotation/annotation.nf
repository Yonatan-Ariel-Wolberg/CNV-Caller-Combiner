nextflow.enable.dsl=2

process DOWNLOAD_ANNOTSV_DB {
    publishDir "${projectDir}/annotsv_db", mode: 'copy'

    output:
    path "AnnotSV_annotations", emit: db_dir

    script:
    def version_tag = params.annotsv_version.startsWith('v') ? params.annotsv_version : "v${params.annotsv_version}"
    def version_clean = params.annotsv_version.replaceAll(/^v/, '')
    def parts = version_clean.split('\\.')
    def major_minor = parts.size() >= 2 ? "${parts[0]}.${parts[1]}" : version_clean
    """
    git clone --branch ${version_tag} --depth 1 https://github.com/lgmgeo/AnnotSV.git
    cd AnnotSV
    
    # Use pre-downloaded human annotations file if available, otherwise download in parallel
    if [ -f "${projectDir}/Annotations_Human_${major_minor}.tar.gz" ]; then
        echo "Using pre-downloaded file from ${projectDir}/Annotations_Human_${major_minor}.tar.gz"
        cp "${projectDir}/Annotations_Human_${major_minor}.tar.gz" .
    else
        echo "Downloading database in parallel..."
        python3 "${projectDir}/bin/parallel_download.py" "https://cstb-icube.fr/~geoffroy/Annotations/Annotations_Human_${major_minor}.tar.gz" "Annotations_Human_${major_minor}.tar.gz" 16
    fi
    
    # Integrity check
    echo "Verifying archive integrity..."
    if ! tar -tzf "Annotations_Human_${major_minor}.tar.gz" > /dev/null; then
        echo "Error: Downloaded file is corrupted." >&2
        exit 1
    fi
    echo "Archive integrity verified."
    
    cp Annotations_Human_${major_minor}.tar.gz Annotations_Human_${version_clean}.tar.gz 2>/dev/null || true
    
    make PREFIX=\$(readlink -f \$(pwd)/../AnnotSV_annotations) install
    make PREFIX=\$(readlink -f \$(pwd)/../AnnotSV_annotations) install-human-annotation
    cd ..
    rm -rf AnnotSV
    """
}

process RUN_ANNOTSV {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(rare_vcf)
    path db_dir

    output:
    tuple val(run_name), path('cohort_rare.annotated.tsv'), emit: annotated_tsv

    script:
    def genome_build = params.get('genome_build', params.get('build', 'GRCh38'))
    """
    export ANNOTSV=\$(readlink -f ${db_dir})
    
    # Run AnnotSV
    AnnotSV \
        -SVinputFile "${rare_vcf}" \
        -outputDir . \
        -outputFile cohort_rare.annotated \
        -genomeBuild ${genome_build} \
        -overwrite 0 || true
        
    # Fallback headers if file not created
    if [ ! -f cohort_rare.annotated.tsv ]; then
        echo -e "Annotation_mode\tSV_chrom\tSV_start\tSV_end\tSV_type\tACMG_class\tAnnotSV_ranking_score\tAnnotSV_ranking_criteria\tGene_name" > cohort_rare.annotated.tsv
    fi
    """
}

process ADD_ACMG_ANNOTATIONS {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_acmg_reports", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), emit: rare_cnvs_annotated
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs_*.csv'), emit: score_separated_csvs

    script:
    """
    python3 "${projectDir}/bin/add-acmg-annotations.py" \
        --input-csv "${rare_cnvs_csv}" \
        --annotsv-tsv "${annotated_tsv}" \
        --output-csv "DDD_Cohort_Rare_CNVs.csv"
    """
}
include { FILTER_RARE_CNVS; FILTER_VCF_BY_RARE_REGIONS; filterRareCnvs; filterVcfByRareRegions } from '../filters/allele-frequency-filters.nf'

process STRATIFY_VCF_BY_ACMG {
    tag "${run_name}"
    label 'pysam'
    publishDir "${params.outdir}/consensus_acmg_vcfs", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_vcf), path(annotated_tsv)

    output:
    tuple val(run_name), path('annotated.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz'), emit: annotated_vcf, optional: true
    tuple val(run_name), path('annotated.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz.tbi'), emit: annotated_tbi, optional: true
    tuple val(run_name), path('ACMG*.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz'), emit: acmg_vcfs, optional: true
    tuple val(run_name), path('ACMG*.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz.tbi'), emit: acmg_tbis, optional: true

    script:
    """
    python3 "${projectDir}/bin/stratify-vcf-by-acmg.py" \\
        --input-vcf "${rare_vcf}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --output-pattern "ACMG{score}.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz"
    """
}

// Backward-compatibility wrappers for legacy pipelines
workflow runAnnotSV {
    take: rare_vcf_ch; db_dir_ch
    main:
        RUN_ANNOTSV(rare_vcf_ch, db_dir_ch)
    emit:
        annotated_tsv = RUN_ANNOTSV.out.annotated_tsv
}

workflow addACMGAnnotations {
    take: acmg_input_ch
    main:
        ADD_ACMG_ANNOTATIONS(acmg_input_ch)
    emit:
        rare_cnvs_annotated = ADD_ACMG_ANNOTATIONS.out.rare_cnvs_annotated
        score_separated_csvs = ADD_ACMG_ANNOTATIONS.out.score_separated_csvs
}

workflow stratifyVcfByAcmg {
    take: vcf_acmg_ch
    main:
        STRATIFY_VCF_BY_ACMG(vcf_acmg_ch)
    emit:
        annotated_vcf = STRATIFY_VCF_BY_ACMG.out.annotated_vcf
        annotated_tbi = STRATIFY_VCF_BY_ACMG.out.annotated_tbi
        acmg_vcfs = STRATIFY_VCF_BY_ACMG.out.acmg_vcfs
        acmg_tbis = STRATIFY_VCF_BY_ACMG.out.acmg_tbis
}
