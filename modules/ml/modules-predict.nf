nextflow.enable.dsl=2

process RUN_CATBOOST_PREDICT {
    label 'train'
    publishDir "${params.outdir}/catboost_prediction", mode: 'copy'

    input:
    path features
    path vcfs
    path truth_labels
    path probes_bed

    output:
    path "predicted_vcfs/*predicted*.vcf.gz", emit: predicted_vcfs, optional: true
    path "predicted_vcfs/*predicted*.vcf.gz.tbi", emit: predicted_tbis, optional: true
    path "predicted_vcfs/catboost_prediction_summary.txt", emit: summary
    path "predicted_vcfs/selected_features.txt", emit: selected_features, optional: true
    path "predicted_vcfs/best_hyperparameters.json", emit: best_hyperparameters, optional: true
    path "predicted_vcfs/catboost_bootstrap_metrics_comparison.md", emit: bootstrap_comparison, optional: true
    path "predicted_vcfs/catboost_bootstrap_metrics.tsv", emit: bootstrap_metrics, optional: true
    path "predicted_vcfs/catboost_shap_values.tsv", emit: shap_values, optional: true
    path "predicted_vcfs/catboost_roc_data.tsv", emit: roc_data, optional: true
    path "predicted_vcfs/catboost_pr_data.tsv", emit: pr_data, optional: true

    script:
    def probes_arg = probes_bed ? "--probes-bed ${probes_bed}" : ""
    def metric_arg = "--eval-metric ${params.get('catboost_eval_metric', 'PRAUC')}"
    def tune_arg = params.get('catboost_tune_hyperparameters', false) ? "--tune-hyperparameters --n-trials ${params.get('catboost_n_trials', 30)}" : ""
    def sel_arg = params.get('catboost_select_features', false) ? "--select-features --n-features ${params.get('catboost_n_features', 12)}" : ""
    """
    mkdir -p predicted_vcfs
    if [ -f "/spaces/ywolberg/nxf_cache/catboost_env/bin/python3" ]; then
        PY_EXEC="/spaces/ywolberg/nxf_cache/catboost_env/bin/python3"
    elif python3 -c "import catboost, pysam, optuna" 2>/dev/null; then
        PY_EXEC="python3"
    else
        echo "Error: CatBoost execution environment not found" >&2
        exit 1
    fi

    \${PY_EXEC} ${projectDir}/bin/train_predict_catboost.py \\
        --features-dir . \\
        --vcf-dir . \\
        --truth-labels ${truth_labels} \\
        --output-dir predicted_vcfs \\
        --threshold ${params.get('catboost_threshold', 0.5)} \\
        ${metric_arg} \\
        ${tune_arg} \\
        ${sel_arg} \\
        ${probes_arg}
    """
}
