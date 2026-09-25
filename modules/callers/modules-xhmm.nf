#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF } from '../common/bgzip-sort-index.nf'

// =====================================================================================
// XHMM MODULE
// =====================================================================================

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params.ref        = ''
params.probes     = ''
params.xhmm_conf  = ''
params.xhmm_batch_size = 50  // BAMs per GATK DepthOfCoverage job; increase for fewer but larger jobs
ref       = file(params.ref       ?: 'NO_FILE', type: 'file')
probes    = file(params.probes    ?: 'NO_FILE', type: 'file')
xhmm_conf = file(params.xhmm_conf ?: 'NO_FILE', type: 'file')
outdir    = file(params.outdir, type: 'dir')

// =====================================================================================
// SPLIT BAM FILES TO GROUPS
// =====================================================================================
// RUN GATK FOR DEPTH OF COVERAGE (FOR INDIVIDUAL SAMPLES):
process GATK_DOC {
    tag { "${sample_id}" }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path("${sample_id}.DATA.sample_interval_summary"), emit: doc_summary

    script:
    mem = Math.max(2, (task.memory.toGiga() * 0.75).toInteger())   // Reserve safe memory headroom for JVM off-heap
    """
    # Resolve symlink to get the real source path of the BAM index
    real_bai=\$(readlink -f $bai 2>/dev/null || echo "$bai")
    rm -f $bai
    cp "\$real_bai" $bai
    touch $bai

    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g -Djava.io.tmpdir=\$PWD" \\
        DepthOfCoverage \\
        -I ${bam} \\
        -L ${probes} \\
        -R ${ref} \\
        --max-depth-per-sample 5000 \\
        --verbosity INFO \\
        --omit-depth-output-at-each-base true \\
        --omit-locus-table true \\
        --min-base-quality 0 \\
        --read-filter MappingQualityReadFilter \\
        --minimum-mapping-quality 20 \\
        --start 1 --stop 5000 --nBins 200 \\
        --include-ref-n-sites true \\
        --count-type COUNT_READS \\
        --output ${sample_id}.DATA
    """
}

// COMBINES GATK DEPTH-OF-COVERAGE OUTPUTS FOR MULTIPLE SAMPLES (AT SAME LOCI):
process COMBINE_DOC {
    tag { 'combine_doc' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(list)
    
    output:
    path("DATA.RD.txt"), emit: combined_doc
    
    script:
    """
    for i in *.sample_interval_summary; do tr ',' '\t' < \$i > \${i%.sample_interval_summary}.fixed_sample_interval_summary ; done
    ls *.fixed_sample_interval_summary > the_list
    xhmm --mergeGATKdepths -o DATA.RD.txt --GATKdepthsList the_list
    """
}

// OPTIONALLY, RUN GATK TO CALCULATE THE PER-TARGET GC CONTENT AND CREATE A LIST OF EXTREME GC CONTENT TARGETS:
process CALC_GC_XHMM {
    tag { 'calc_gc' }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    output:
    path("DATA.locus_GC.txt"), emit: gc_targets
    path("extreme_gc_targets.txt"), emit: extreme_gc_targets
    
    script:
    """
    gatk AnnotateIntervals \
        -R ${ref} \
        -L ${probes} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O gc.tmp

    grep -v "^@" gc.tmp | grep -v "^CONTIG" | awk '{ print \$1":"\$2"-"\$3"\t"\$4 }' > DATA.locus_GC.txt
    cat DATA.locus_GC.txt | awk '{ if (\$2 < 0.1 || \$2 > 0.9) print \$1 }' > extreme_gc_targets.txt
    """
}

// FILTERS SAMPLES AND TARGETS AND THEN MEAN-CENTERS THE TARGETS:
process FILTER_SAMPLES {
    tag { 'filter_samples' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true    

    input:
    path(combined_doc)
    path(extreme_gc_targets)
    
    output:
    path("DATA.filtered_centered.RD.txt"), emit: filtered_centered
    path("DATA.filtered_centered.RD.txt.filtered_targets.txt"), emit: excluded_filtered_targets
    path("DATA.filtered_centered.RD.txt.filtered_samples.txt"), emit: excluded_filtered_samples
    
    script:
    """
    xhmm --matrix -r DATA.RD.txt --centerData --centerType target \
        -o DATA.filtered_centered.RD.txt \
        --outputExcludedTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeTargets extreme_gc_targets.txt \
        --minTargetSize 10 --maxTargetSize 10000 \
        --minMeanTargetRD 10 --maxMeanTargetRD 500 \
        --minMeanSampleRD 25 --maxMeanSampleRD 200 \
        --maxSdSampleRD 150
   """
}

// RUNS PCA ON MEAN-CENTERED DATA:
process RUN_PCA {
    tag { 'run_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    
    output:
    tuple val("pca_data"), path("DATA.RD_PCA*"), emit: pca_data

    script:
    """
    xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA
    """
}

// NORMALIZES MEAN-CENTERED DATA USING PCA INFORMATION:
process NORMALISE_PCA {
    tag { 'norm_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    tuple val(pca), path(pca_data)
    
    output:
    path("DATA.PCA_normalized.txt"), emit: data_pca_norm

    script:
    """
    xhmm --normalize -r DATA.filtered_centered.RD.txt \
        --PCAfiles DATA.RD_PCA \
        --normalizeOutput DATA.PCA_normalized.txt \
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
    """
}

// FILTERS AND Z-SCORE CENTERS (BY SAMPLE) THE PCA-NORMALIZED DATA:
process FILTER_ZSCORE {
    tag { 'filter_zscore' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(data_pca_norm)
    
    output:
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt"), emit: pca_norm_zscore
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"), emit: excluded_zscore_targets
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"), emit: excluded_zscore_samples
    
    script:
    """
    xhmm --matrix -r DATA.PCA_normalized.txt \
        --centerData --centerType sample --zScoreData \
        -o DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        --outputExcludedTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        --maxSdTargetRD 30
    """
}

// FILTERS ORIGINAL READ-DEPTH DATA TO BE THE SAME AS FILTERED, NORMALIZED DATA:
process FILTER_RD {
    tag { 'filter_rd' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(combined_doc)
    path(excluded_filtered_targets)
    path(excluded_zscore_targets)
    path(excluded_filtered_samples)
    path(excluded_zscore_samples)
    
    output:
    path("DATA.same_filtered.RD.txt"), emit: orig_filtered

    script:
    """
    xhmm --matrix -r DATA.RD.txt \
        --excludeTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --excludeTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --excludeSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        -o DATA.same_filtered.RD.txt
    """
}

// DISCOVERS CNVS IN NORMALIZED DATA:
process DISCOVER_CNVS {
    tag { 'discover_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    
    output:
    path("*"), emit: cnvs
    
    script:
    """
    xhmm --discover -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -c DATA.xcnv -a DATA.aux_xcnv -s DATA
    """
}

// GENOTYPES DISCOVERED CNVS IN ALL SAMPLES:
process GENOTYPE_CNVS {
    tag { 'genotype_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    path(cnvs)
    
    output:
    path("*"), emit: genotypes
    
    script:
    """
    xhmm --genotype -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -g DATA.xcnv -F ${ref} \
        -v DATA.vcf
    """
}

// SPLITS COMBINED VCF INTO INDIVIDUAL SAMPLE VCFs
process SPLIT_VCF {
    tag { 'split_vcf' }
    label 'bcftools'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path vcf_file
    path fasta_fai
    
    output:
    path("*_XHMM_output.vcf"), emit: individual_vcfs
    
    script:
    """
    # Generate contig headers from the FAI file
    awk '{printf("##contig=<ID=%s,length=%s>\\n", \$1, \$2)}' ${fasta_fai} > contigs.hdr

    # Inject the contig headers right after the ##fileformat line in the VCF
    sed -e '/##fileformat/r contigs.hdr' ${vcf_file} > temp.vcf

    # Split the multi-sample VCF by sample using bcftools view -s and awk.
    # Keeping only variants actually called for each sample (GT != ref and GT != missing),
    # and ensure the CHROM column carries the 'chr' prefix required by downstream tools.
    for sample in \$(bcftools query -l temp.vcf); do
        sample_clean=\$(echo "\${sample}" | sed 's/\\.md\$//')
        bcftools view -s "\${sample}" -i 'GT!="." && GT!="0"' temp.vcf | \\
        awk -v smp_clean="\${sample_clean}" '
        BEGIN { OFS="\\t" }
        /^##/ { print; next }
        /^#CHROM/ {
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, smp_clean
            next
        }
        {
            chrom = \$1
            if (chrom !~ /^chr/) {
                chrom = "chr" chrom
            }
            split(\$10, a, ":")
            gt = a[1]
            alt_allele = ""
            if (gt == "1") {
                alt_allele = "<DEL>"
            } else if (gt == "2") {
                alt_allele = "<DUP>"
            }
            
            if (alt_allele != "") {
                sub(/^[^:]+/, gt, \$10)
                print chrom, \$2, \$3, \$4, alt_allele, \$6, \$7, \$8, \$9, \$10
            } else {
                print chrom, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10
            }
        }
        ' > \${sample_clean}_XHMM_output.vcf
    done
    rm -f temp.vcf contigs.hdr
    """
}


// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow XHMM {
    take:
    bam_ch  // channel: tuple(sample_id, bam, bai)

    main:
    // Step 1: Calculate depth of coverage per sample in parallel
    GATK_DOC(bam_ch)

    // Step 2: Combine all individual depth-of-coverage outputs
    COMBINE_DOC(GATK_DOC.out.doc_summary.collect())

    // Step 4: Calculate per-target GC content and identify extreme GC targets
    CALC_GC_XHMM()

    // Step 5: Filter samples and targets, then mean-center the data
    FILTER_SAMPLES(COMBINE_DOC.out.combined_doc, CALC_GC_XHMM.out.extreme_gc_targets)

    // Step 6: Run PCA on mean-centered data
    RUN_PCA(FILTER_SAMPLES.out.filtered_centered)

    // Step 7: Normalize mean-centered data using PCA information
    NORMALISE_PCA(FILTER_SAMPLES.out.filtered_centered, RUN_PCA.out.pca_data)

    // Step 8: Filter and z-score center the PCA-normalized data
    FILTER_ZSCORE(NORMALISE_PCA.out.data_pca_norm)

    // Step 9: Filter original read-depth data to match the normalized filtered data
    FILTER_RD(
        COMBINE_DOC.out.combined_doc,
        FILTER_SAMPLES.out.excluded_filtered_targets,
        FILTER_ZSCORE.out.excluded_zscore_targets,
        FILTER_SAMPLES.out.excluded_filtered_samples,
        FILTER_ZSCORE.out.excluded_zscore_samples
    )

    // Step 10: Discover CNVs in the normalized data
    DISCOVER_CNVS(FILTER_RD.out.orig_filtered, FILTER_ZSCORE.out.pca_norm_zscore)

    // Step 11: Genotype discovered CNVs across all samples
    GENOTYPE_CNVS(
        FILTER_RD.out.orig_filtered,
        FILTER_ZSCORE.out.pca_norm_zscore,
        DISCOVER_CNVS.out.cnvs
    )

    // Step 12: Split the combined VCF into per-sample VCFs
    vcf_ch = GENOTYPE_CNVS.out.genotypes
        .flatten()
        .filter { it.name == 'DATA.vcf' }

    SPLIT_VCF(vcf_ch, file((params.ref ?: 'NO_FILE') + '.fai'))

    // Step 13: Compress, sort, index, and annotate raw per-sample VCFs with TOOL=XHMM
    // (Quality score filtering is handled centrally in quality-filters.nf)
    BGZIP_SORT_INDEX_VCF(SPLIT_VCF.out.individual_vcfs.flatten(), 'XHMM', 'out_XHMM')

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
