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
    path "predicted_vcfs/*predicted*.vcf*", emit: predicted_vcfs, optional: true
    path "predicted_vcfs/catboost_prediction_summary.txt", emit: summary
    path "predicted_vcfs/catboost_bootstrap_metrics_comparison.md", emit: bootstrap_comparison, optional: true
    path "predicted_vcfs/catboost_bootstrap_metrics.tsv", emit: bootstrap_metrics, optional: true
    path "predicted_vcfs/catboost_shap_values.tsv", emit: shap_values, optional: true
    path "predicted_vcfs/catboost_roc_data.tsv", emit: roc_data, optional: true
    path "predicted_vcfs/catboost_pr_data.tsv", emit: pr_data, optional: true

    script:
    def probes_arg = probes_bed ? "--probes-bed ${probes_bed}" : ""
    """
    mkdir -p predicted_vcfs
    python3 -c "import catboost, pysam" 2>/dev/null || mamba install -y -n base -c conda-forge catboost xgboost scikit-learn imbalanced-learn pandas numpy pysam

    python3 ${projectDir}/bin/train_predict_catboost.py \
        --features-dir . \
        --vcf-dir . \
        --truth-labels ${truth_labels} \
        --output-dir predicted_vcfs \
        --threshold ${params.get('catboost_threshold', 0.5)} \
        ${probes_arg}
    """
}
