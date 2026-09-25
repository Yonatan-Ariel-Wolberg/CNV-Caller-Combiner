#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import all process modules
include { CANOES } from '../callers/modules-canoes.nf'
include { XHMM } from '../callers/modules-xhmm.nf'
include { CLAMMS } from '../callers/modules-clamms.nf'
include { DRAGEN } from '../callers/modules-icav2-dragen.nf'
include { CNVKIT } from '../callers/modules-cnvkit.nf'
include { GATK_GCNV } from '../callers/modules-gatk-gcnv.nf'
include { EXOMEDEPTH } from '../callers/modules-exomedepth.nf'
include { SURVIVOR } from '../sv-mergers/modules-survivor.nf'
include { TRUVARI } from '../sv-mergers/modules-truvari.nf'
include { FEATURE_EXTRACTION } from '../ml/modules-feature-extraction.nf'
include { NORMALISE } from '../normalise/modules-normalise.nf'
include { TRAIN } from '../ml/modules-train.nf'
include { EVALUATE; EVALUATE_COMPARATIVE; EVALUATE_ACMG_CSV; PLOT_EVALUATION_METRICS } from '../evaluation/performance-metric-evaluation.nf'
include { SAMTOOLS_QC_BAM; SAMTOOLS_QC_CRAM; QUALIMAP; MULTIQC; INFER_SEX_AND_PLOIDY; ESTIMATE_CONTAMINATION } from '../sample-qc/target-aware-exome-sequencing-quality-assess.nf'
include { SVABA; RUN_SVABA } from '../callers/modules-svaba.nf'
include { REFINE_CONSENSUS_BREAKPOINTS; RUN_BREAKPOINT_REFINEMENT } from '../refinement/breakpoint-refinement.nf'

include { MERGE_CONSENSUS_COHORT; SPLIT_COHORT_VCF; SANITIZE_COHORT_VCF; BCFTOOLS_FLATTEN_COHORT; sanitizeCohortVcf; bcftoolsFlattenCohort } from '../consensus/consensus-postprocess.nf'
include { DOWNLOAD_ANNOTSV_DB; RUN_ANNOTSV; ADD_ACMG_ANNOTATIONS; STRATIFY_VCF_BY_ACMG; runAnnotSV; addACMGAnnotations; stratifyVcfByAcmg } from '../annotation/annotation.nf'
include { FILTER_RARE_CNVS; FILTER_VCF_BY_RARE_REGIONS; filterRareCnvs; filterVcfByRareRegions } from '../filters/allele-frequency-filters.nf'
include { ANALYZE_CNV_QUALITY_METRICS } from '../evaluation/quality-metric-evaluation.nf'
include { RUN_CATBOOST_PREDICT } from '../ml/modules-predict.nf'
include { CONVERT_CRAM_TO_BAM_VERIFIED; CONVERT_BAM_TO_CRAM_VERIFIED; RUN_CRAM_TO_BAM_CONVERSION; RUN_BAM_TO_CRAM_CONVERSION; CONVERT_CRAM_TO_TARGET_PADDED_BAM; RUN_CRAM_TO_TARGET_PADDED_BAM_CONVERSION } from '../conversions/modules-conversions.nf'

workflow RUN_CANOES {
    take: 
        bams
        chroms
        fai
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        CANOES(bam_list, fai, Channel.from(chroms))
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CANOES.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CANOES'
            )
        }
}

workflow RUN_XHMM {
    take: 
        bams
    main:
        XHMM(bams)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                XHMM.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'XHMM'
            )
        }
}

workflow RUN_CLAMMS {
    take: 
        bams
        fai
    main:
        sample_list = bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
        CLAMMS(bams, fai, sample_list)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CLAMMS.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CLAMMS'
            )
        }
}

workflow RUN_DRAGEN {
    take: 
        cramPairs
    main:
        DRAGEN(cramPairs)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                DRAGEN.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'DRAGEN'
            )
        }
}

workflow RUN_CNVKIT {
    take: 
        bams
        fasta
        targets
        refflat
    main:
        // Guard against empty BAM input so workflow fails fast with a clear error.
        bams_nonempty = bams.ifEmpty { error "RUN_CNVKIT received no BAMs. Check --bams input/glob." }
        CNVKIT(bams_nonempty, fasta, targets, refflat, bams_nonempty.first().map { it[1] })
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CNVKIT.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CNVKIT'
            )
        }
}

workflow RUN_GCNV {
    take: 
        bams
        fasta
        fai
        dict
        targets
    main:
        GATK_GCNV(bams, fasta, fai, dict, targets)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                GATK_GCNV.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'GCNV'
            )
        }
}

workflow RUN_EXOMEDEPTH {
    take:
        bams
        fai
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        probes_f = file(params.get('probes', params.get('targets', params.get('exome_targets', params.get('probes_bed', '')))))
        EXOMEDEPTH(bam_list, probes_f, fai)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                EXOMEDEPTH.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'EXOMEDEPTH'
            )
        }
}

workflow RUN_EVALUATE {
    take:
        vcf_ch
        truth_bed
        probes_bed
        caller_name
    main:
        EVALUATE(vcf_ch, truth_bed, probes_bed, caller_name)
}

workflow RUN_SURVIVOR {
    take: 
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        SURVIVOR(grouped_vcfs)
    emit:
        union_vcf = SURVIVOR.out.union_vcf
}

workflow RUN_TRUVARI {
    take: 
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        TRUVARI(grouped_vcfs)
    emit:
        merged_vcf = TRUVARI.out.merged_vcf
}

workflow RUN_FEATURE_EXTRACTION {
    take:
        // Channel of tuples:
        // [ sample_id, merged_vcf, collapsed_vcf, tool_vcfs_str, merger_mode,
        //   bam_file, reference_fasta, bed_file, mappability_file ]
        // collapsed_vcf should be [] when not available (SURVIVOR mode or absent).
        feature_inputs_ch
    main:
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_NORMALISE {
    take:
        // Channel of tuples: [ sample_name, vcf_path, caller_name ]
        vcf_ch
    main:
        NORMALISE(vcf_ch)
}

workflow RUN_TRAIN {
    take:
        // Channel of individual *_features.tsv file paths (one per sample).
        features_tsv_ch
        // Single-element channel containing the path to the truth-labels TSV.
        // Required columns: sample_id, chrom, start, end, cnv_type, truth_label
        truth_labels_ch
        // Optional single-element channel containing probes BED.
        // Used for probe-overlap fallback when feature/truth CNV coordinates differ.
        probes_bed_ch
    main:
        TRAIN(features_tsv_ch, truth_labels_ch, probes_bed_ch)
}

workflow RUN_SURVIVOR_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        SURVIVOR(grouped_vcfs)
        feature_inputs_ch = SURVIVOR.out.union_vcf
            .map { sample_id, merged_vcf -> build_feature_inputs(sample_id, merged_vcf, []) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_TRUVARI_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        TRUVARI(grouped_vcfs)
        feature_inputs_ch = TRUVARI.out.merged_vcf
            .join(TRUVARI.out.collapsed_vcf)
            .map { sample_id, merged_vcf, collapsed_vcf -> build_feature_inputs(sample_id, merged_vcf, collapsed_vcf) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_QC {
    take:
        input_dir
        output_dir
        ref_genome
        target_bed
        mapq_threshold
        threads
    main:
        input_ch = Channel.fromPath("${input_dir}/*.{bam,cram}")

        input_ch.branch {
            bams: it.name.endsWith('.bam')
            crams: it.name.endsWith('.cram')
        }.set { branched_files }

        // Channel for optional target bed file
        bed_ch = target_bed ? Channel.fromPath(target_bed).collect() : Channel.value([])

        // Execute BAM pipeline
        SAMTOOLS_QC_BAM(branched_files.bams, bed_ch)
        
        // Execute CRAM pipeline
        if (ref_genome) {
            ref_ch = Channel.fromPath(ref_genome).collect()
            SAMTOOLS_QC_CRAM(branched_files.crams, ref_ch, bed_ch)
        }

        // QualiMap
        bam_pairs = SAMTOOLS_QC_BAM.out.bam_and_index
        cram_pairs = ref_genome ? SAMTOOLS_QC_CRAM.out.cram_and_index : Channel.empty()
        qualimap_input = bam_pairs.mix(cram_pairs)

        QUALIMAP(qualimap_input, bed_ch)

        // Collect all raw reports for MultiQC
        bam_metrics = SAMTOOLS_QC_BAM.out.flagstat.mix(SAMTOOLS_QC_BAM.out.stats).mix(SAMTOOLS_QC_BAM.out.idxstats)
        cram_metrics = ref_genome ? SAMTOOLS_QC_CRAM.out.flagstat.mix(SAMTOOLS_QC_CRAM.out.stats).mix(SAMTOOLS_QC_CRAM.out.idxstats) : Channel.empty()

        all_qc_metrics = bam_metrics.mix(cram_metrics).mix(QUALIMAP.out.qc_dir).collect()

        MULTIQC(all_qc_metrics)

        // Automated Sex & Ploidy Inference
        idxstats_ch = SAMTOOLS_QC_BAM.out.idxstats.mix(ref_genome ? SAMTOOLS_QC_CRAM.out.idxstats : Channel.empty()).collect()
        INFER_SEX_AND_PLOIDY(idxstats_ch)

        // Automated Reference-Free Contamination Estimator (Method A)
        ESTIMATE_CONTAMINATION(branched_files.bams.collect(), bed_ch)
}

workflow RUN_CONSENSUS_ANNOTATION {
    take:
        run_name
        sample_vcfs_ch
        annotsv_db_ch
        probes_bed
    main:
        // Cohort-wide merge
        MERGE_CONSENSUS_COHORT(sample_vcfs_ch.collect())
        
        cohort_ch = Channel.of(run_name).combine(MERGE_CONSENSUS_COHORT.out.merged_vcf)
        
        bcftoolsFlattenCohort(cohort_ch)
        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
        
        slice_joined = cohort_ch.join(filterRareCnvs.out.rare_data)
        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
        
        runAnnotSV(filterVcfByRareRegions.out.rare_vcf, annotsv_db_ch)
        
        acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })

        vcf_acmg_joined = filterVcfByRareRegions.out.rare_vcf.join(runAnnotSV.out.annotated_tsv)
        stratifyVcfByAcmg(vcf_acmg_joined)

        // Quality assessment
        metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
            .join(runAnnotSV.out.annotated_tsv)
            .combine(sample_vcfs_ch.collect().map { [ it ] })
            .map { run_id, rare_csv, annot_tsv, vcfs ->
                [ run_id, rare_csv, annot_tsv, vcfs ]
            }
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_inputs,
            file(probes_bed)
        )

        // ACMG Evaluators
        if (params.get('map_file', false) && params.get('truth_bed', false)) {
            eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
            eval_joined = eval_inputs.combine(Channel.of(probes_bed))
            
            EVALUATE_ACMG_CSV(
                eval_joined.map { [ it[0], it[1], it[2] ] },
                file(params.map_file),
                file(params.truth_bed)
            )

            PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
        }
}
