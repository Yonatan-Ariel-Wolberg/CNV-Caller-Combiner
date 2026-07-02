#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESS: TRAIN XGBOOST CLASSIFIER
// =====================================================================================
//
// Trains a gradient-boosted tree classifier (XGBoost) on the CNV feature
// matrices produced by the feature_extraction workflow.
//
// All *_features.tsv files staged into the process working directory are
// discovered automatically, so the caller simply collects every
// EXTRACT_FEATURES output into a single channel before calling this process.
//
// Input tuple layout:
//   path  features_tsv_files -- collected list of *_features.tsv files
//   path  truth_labels       -- TSV with columns: sample_id, chrom, start, end,
//                               cnv_type, truth_label  (1 = true CNV, 0 = false positive)
//   val   probes_bed         -- optional BED file path (CHR, START, END) used for
//                               probe-overlap fallback matching when CNV
//                               coordinates differ between features/truth labels
//
// Outputs:
//   model  -- trained XGBoost model saved in JSON format (cnv_model.json)
//   report -- plain-text training summary (training_report.txt)
//   roc_plot / pr_plot -- ROC and precision-recall curve plots (SVG)
//   roc_data / pr_data -- ROC and precision-recall points (TSV)
//   shap_values -- per-call SHAP values (TSV)
//   shap_summary_plot / shap_beeswarm_plot -- SHAP summary plots (SVG)

process TRAIN_XGBOOST {
    tag "train"
    label 'train'
    publishDir "${outdir}/out_TRAIN", mode: 'copy', overwrite: true

    input:
    path(features_tsv_files)
    path(truth_labels)
    val(probes_bed)

    output:
    path("cnv_model.json"),      emit: model
    path("training_report.txt"), emit: report
    path("roc_curve.svg"),       emit: roc_plot
    path("pr_curve.svg"),        emit: pr_plot
    path("roc_curve.tsv"),       emit: roc_data
    path("pr_curve.tsv"),        emit: pr_data
    path("shap_values.tsv"),     emit: shap_values
    path("shap_summary_bar.svg"), emit: shap_summary_bar_plot
    path("shap_summary_beeswarm.svg"), emit: shap_beeswarm_plot

    script:
    // All feature TSV files are staged into the working directory by Nextflow.
    // The Python script globs *_features.tsv from the current directory ('.').
    def probes_escaped = probes_bed ? probes_bed.toString().replace("'", "'\\''") : null
    def probes_arg = probes_escaped ? "--probes_bed '${probes_escaped}'" : ""
    """
    set -euo pipefail

    mamba install -y -n base -c conda-forge \
        xgboost=2.1.4 \
        scikit-learn=1.6.1 \
        imbalanced-learn=0.13.0 \
        pandas=2.2.3 \
        numpy=2.2.3

    python ${projectDir}/bin/train_xgboost.py \
        --features_dir   '.' \
        --truth_labels   '${truth_labels}' \
        ${probes_arg} \
        --output_model   cnv_model.json \
        --output_report  training_report.txt \
        --output_roc_plot roc_curve.svg \
        --output_pr_plot  pr_curve.svg \
        --output_roc_data roc_curve.tsv \
        --output_pr_data  pr_curve.tsv \
        --output_shap_values shap_values.tsv \
        --output_shap_summary_plot shap_summary_bar.svg \
        --output_shap_beeswarm_plot shap_summary_beeswarm.svg
    """
}

// =====================================================================================
// SUB-WORKFLOW
// =====================================================================================

workflow TRAIN {
    take:
        // Channel of individual *_features.tsv file paths (one per sample).
        // Internally collected before passing to TRAIN_XGBOOST so that the
        // classifier sees all samples in a single training run.
        features_tsv_ch
        // Channel containing a single path: the truth-labels TSV.
        truth_labels_ch
        // Optional channel containing a single probes BED path.
        probes_bed_ch

    main:
        features_tsv_ch.ifEmpty { error "TRAIN received no feature TSV files. Check that features_tsv_ch is not empty." }
        TRAIN_XGBOOST(features_tsv_ch.collect(), truth_labels_ch, probes_bed_ch)

    emit:
        model              = TRAIN_XGBOOST.out.model
        report             = TRAIN_XGBOOST.out.report
        roc_plot           = TRAIN_XGBOOST.out.roc_plot
        pr_plot            = TRAIN_XGBOOST.out.pr_plot
        roc_data           = TRAIN_XGBOOST.out.roc_data
        pr_data            = TRAIN_XGBOOST.out.pr_data
        shap_values        = TRAIN_XGBOOST.out.shap_values
        shap_summary_bar_plot = TRAIN_XGBOOST.out.shap_summary_bar_plot
        shap_beeswarm_plot    = TRAIN_XGBOOST.out.shap_beeswarm_plot
}
