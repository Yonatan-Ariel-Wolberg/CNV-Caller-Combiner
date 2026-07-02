#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: train_compare_models.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
"""
train_compare_models.py

Trains and compares 5 machine learning classification models
(CatBoost vs LightGBM vs XGBoost vs HistGradientBoosting vs Random Forest)
to distinguish true-positive CNV calls from false-positive background noise.

Implements GroupKFold cross-validation (grouped by sample_id to prevent leakage)
and bootstrap resampling (B >= 1000) for 95% confidence intervals of sensitivity,
precision, specificity, F2, F0.5, and MCC. Outputs comparative curves and SHAP
beeswarm plots.
"""

import os
import sys
import glob
import math
import argparse
import logging
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.model_selection import GroupKFold
from sklearn.metrics import precision_recall_curve, roc_curve, auc, fbeta_score, matthews_corrcoef, confusion_matrix
from sklearn.impute import SimpleImputer
from sklearn.utils import resample
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier

# Ignore warnings
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="Compare 5 ML models for CNV classification.")
    parser.add_argument('--features-dir', required=True,
                        help="Directory containing per-sample feature TSVs (from EXTRACT_FEATURES).")
    parser.add_argument('--truth-bed', required=True,
                        help="Path to Sanger truthset BED file.")
    parser.add_argument('--map-file', required=True,
                        help="Path to EGAF-EGAN map file.")
    parser.add_argument('--output-dir', required=True,
                        help="Directory to save plots and reports.")
    parser.add_argument('--n-bootstraps', type=int, default=1000,
                        help="Number of bootstrap resamples for confidence intervals.")
    return parser.parse_args()

def load_map_file(map_file):
    """Load EGAF to EGAN Sample ID mapping."""
    egaf_to_egan = {}
    with open(map_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                egaf_to_egan[parts[0]] = parts[1]
    return egaf_to_egan

def load_truth_bed(truth_file):
    """Load 5-column truth BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)."""
    df = pd.read_csv(truth_file, sep='\t', header=None, comment='#', keep_default_na=False)
    if df.shape[1] < 5:
        for col_idx in range(df.shape[1], 5):
            df[col_idx] = ""
    df = df.iloc[:, [0, 1, 2, 3, 4]]
    df.columns = ['chr', 'start', 'end', 'cnv_type', 'sample']
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    df = df.dropna(subset=['start', 'end'])
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['chr'] = df['chr'].astype(str)
    df['sample'] = df['sample'].astype(str).str.strip()
    return df

def map_calls_to_truth(calls_df, truth_df, egaf_to_egan):
    """
    Map candidate calls to the Sanger truthset.
    Assigns class label 1 if the call overlaps with a truth interval
    on the same sample and chromosome, else 0.
    """
    labels = []
    
    # Pre-group truth by sample and chromosome for O(1) lookups
    truth_lookup = {}
    for (sample, chrom), grp in truth_df.groupby(['sample', 'chr']):
        truth_lookup[(sample, chrom)] = list(zip(grp['start'].values, grp['end'].values))

    for idx, row in calls_df.iterrows():
        sample_egaf = str(row['sample_id']).strip()
        sample_egan = egaf_to_egan.get(sample_egaf, sample_egaf)
        chrom = str(row['chrom']).strip()
        # Ensure 'chr' prefix matches
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
            
        start = int(row['start'])
        end = int(row['end'])
        
        # Check overlaps
        is_tp = 0
        truth_intervals = truth_lookup.get((sample_egan, chrom), [])
        for t_start, t_end in truth_intervals:
            if t_start < end and t_end > start:
                is_tp = 1
                break
        labels.append(is_tp)
        
    return labels

def compute_bootstrap_metrics(y_true, y_pred, y_prob, n_bootstraps=1000):
    """
    Compute Sensitivity, Specificity, Precision, F2, F0.5, and MCC
    with 95% Confidence Intervals using bootstrap resampling.
    """
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    y_prob = np.array(y_prob)
    
    metrics = {m: [] for m in ['sens', 'spec', 'prec', 'f2', 'f0.5', 'mcc', 'auc_roc', 'auc_pr']}
    
    for _ in range(n_bootstraps):
        idx = np.random.choice(len(y_true), len(y_true), replace=True)
        y_t_b, y_pred_b, y_prob_b = y_true[idx], y_pred[idx], y_prob[idx]
        
        # Confusion matrix
        cm = confusion_matrix(y_t_b, y_pred_b, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f2 = fbeta_score(y_t_b, y_pred_b, beta=2, zero_division=0)
        f0_5 = fbeta_score(y_t_b, y_pred_b, beta=0.5, zero_division=0)
        mcc = matthews_corrcoef(y_t_b, y_pred_b)
        
        # ROC AUC
        fpr, tpr, _ = roc_curve(y_t_b, y_prob_b)
        auc_roc = auc(fpr, tpr) if len(np.unique(y_t_b)) > 1 else 0.5
        
        # PR AUC
        p, r, _ = precision_recall_curve(y_t_b, y_prob_b)
        auc_pr = auc(r, p) if len(np.unique(y_t_b)) > 1 else 0.0
        
        metrics['sens'].append(sens)
        metrics['spec'].append(spec)
        metrics['prec'].append(prec)
        metrics['f2'].append(f2)
        metrics['f0.5'].append(f0_5)
        metrics['mcc'].append(mcc)
        metrics['auc_roc'].append(auc_roc)
        metrics['auc_pr'].append(auc_pr)
        
    summary = {}
    for k, v in metrics.items():
        summary[k] = {
            'mean': float(np.mean(v)),
            'ci_lower': float(np.percentile(v, 2.5)),
            'ci_upper': float(np.percentile(v, 97.5))
        }
    return summary

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. Load mappings and truthset
    logging.info("Loading sample mapping and truthset...")
    egaf_to_egan = load_map_file(args.map_file)
    truth_df = load_truth_bed(args.truth_bed)
    
    # 2. Load all feature TSVs
    feature_files = glob.glob(os.path.join(args.features_dir, "**/*_features.tsv"), recursive=True)
    if not feature_files:
        logging.error(f"No feature files found in {args.features_dir}!")
        sys.exit(1)
        
    logging.info(f"Found {len(feature_files)} feature files. Concatenating and mapping truth labels...")
    dfs = []
    for f in feature_files:
        temp_df = pd.read_csv(f, sep='\t')
        if not temp_df.empty:
            dfs.append(temp_df)
            
    if not dfs:
        logging.error("All feature files are empty!")
        sys.exit(1)
        
    master_df = pd.concat(dfs, ignore_index=True)
    
    # Extract ground truth labels
    master_df['label'] = map_calls_to_truth(master_df, truth_df, egaf_to_egan)
    
    # Report class balance
    class_counts = master_df['label'].value_counts()
    pos_count = class_counts.get(1, 0)
    neg_count = class_counts.get(0, 0)
    logging.info(f"Dataset summary: {len(master_df)} calls. Class 1 (TP)={pos_count}, Class 0 (FP)={neg_count}")
    
    if pos_count == 0:
        logging.error("No True Positive calls found in the dataset! Model training cannot proceed.")
        sys.exit(1)
        
    # Calculate scale pos weight for imbalanced classifiers
    scale_pos_weight = float(neg_count) / float(pos_count) if pos_count > 0 else 1.0
    
    # Define features to exclude (metadata / targets)
    exclude_cols = {'sample_id', 'chrom', 'start', 'end', 'label', 'chrom_encoded'}
    features_list = [c for c in master_df.columns if c not in exclude_cols]
    logging.info(f"Training on {len(features_list)} features: {features_list}")
    
    X = master_df[features_list].astype(np.float32)
    y = master_df['label'].values
    groups = master_df['sample_id'].values
    
    # RF and some methods don't support NaNs. We impute with mean for RF
    imputer = SimpleImputer(strategy='mean')
    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=features_list)
    
    # 3. Model definitions
    classifiers = {
        'Random Forest': RandomForestClassifier(n_estimators=200, max_depth=10, class_weight='balanced', random_state=42),
        'HistGradientBoosting': HistGradientBoostingClassifier(max_iter=150, max_leaf_nodes=31, class_weight='balanced', random_state=42),
        'XGBoost': XGBClassifier(n_estimators=150, max_depth=6, scale_pos_weight=scale_pos_weight, random_state=42),
        'LightGBM': LGBMClassifier(n_estimators=150, max_depth=6, scale_pos_weight=scale_pos_weight, random_state=42, verbose=-1),
        'CatBoost': CatBoostClassifier(iterations=150, depth=6, scale_pos_weight=scale_pos_weight, random_state=42, verbose=0)
    }
    
    # Initialize out-of-fold predictions
    oof_probs = {model_name: np.zeros(len(y)) for model_name in classifiers}
    oof_preds = {model_name: np.zeros(len(y)) for model_name in classifiers}
    
    # 4. GroupKFold Cross-Validation
    logging.info("Running 5-fold GroupKFold Cross Validation...")
    gkf = GroupKFold(n_splits=5)
    
    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups=groups)):
        logging.info(f"Fold {fold + 1}...")
        
        # Imputed splits for RF
        X_train_imp, X_val_imp = X_imputed.iloc[train_idx], X_imputed.iloc[val_idx]
        # Native NaN splits for Boosting models
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        for name, clf in classifiers.items():
            if name == 'Random Forest':
                clf.fit(X_train_imp, y_train)
                oof_probs[name][val_idx] = clf.predict_proba(X_val_imp)[:, 1]
                oof_preds[name][val_idx] = clf.predict(X_val_imp)
            else:
                clf.fit(X_train, y_train)
                oof_probs[name][val_idx] = clf.predict_proba(X_val)[:, 1]
                oof_preds[name][val_idx] = clf.predict(X_val)
                
    # 5. Evaluate and compute bootstrap metrics
    logging.info("Computing validation metrics with 95% Confidence Intervals...")
    model_metrics = {}
    for name in classifiers:
        model_metrics[name] = compute_bootstrap_metrics(y, oof_preds[name], oof_probs[name], args.n_bootstraps)
        
    # Baseline: concordance >= 2 (if concordance feature exists)
    if 'concordance' in features_list:
        logging.info("Evaluating baseline model (concordance >= 2)...")
        baseline_pred = (X['concordance'] >= 2.0).astype(int).values
        baseline_prob = X['concordance'].values / float(max(1.0, X['concordance'].max()))
        model_metrics['Baseline (Concordance >= 2)'] = compute_bootstrap_metrics(y, baseline_pred, baseline_prob, args.n_bootstraps)
        
    # 6. Generate comparative plots
    logging.info("Generating comparative ROC and Precision-Recall plots...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Colors
    colors = {
        'Random Forest': '#4a86e8',
        'HistGradientBoosting': '#8e7cc3',
        'XGBoost': '#e06666',
        'LightGBM': '#ffd966',
        'CatBoost': '#93c47d',
        'Baseline (Concordance >= 2)': '#999999'
    }
    
    # Plot ROC curves
    ax_roc = axes[0]
    for name in model_metrics:
        if name in colors:
            prob = oof_probs[name] if name != 'Baseline (Concordance >= 2)' else baseline_prob
            fpr, tpr, _ = roc_curve(y, prob)
            roc_auc = auc(fpr, tpr)
            ax_roc.plot(fpr, tpr, label=f"{name} (AUC = {roc_auc:.3f})", color=colors[name], lw=2)
    ax_roc.plot([0, 1], [0, 1], 'k--', lw=1.5)
    ax_roc.set_xlim([0.0, 1.0])
    ax_roc.set_ylim([0.0, 1.05])
    ax_roc.set_xlabel('False Positive Rate')
    ax_roc.set_ylabel('True Positive Rate')
    ax_roc.set_title('ROC Curves')
    ax_roc.legend(loc="lower right")
    
    # Plot Precision-Recall curves
    ax_pr = axes[1]
    for name in model_metrics:
        if name in colors:
            prob = oof_probs[name] if name != 'Baseline (Concordance >= 2)' else baseline_prob
            p, r, _ = precision_recall_curve(y, prob)
            pr_auc = auc(r, p)
            ax_pr.plot(r, p, label=f"{name} (AUC = {pr_auc:.3f})", color=colors[name], lw=2)
    ax_pr.set_xlim([0.0, 1.0])
    ax_pr.set_ylim([0.0, 1.05])
    ax_pr.set_xlabel('Recall')
    ax_pr.set_ylabel('Precision')
    ax_pr.set_title('Precision-Recall Curves')
    ax_pr.legend(loc="lower left")
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'comparative_performance_curves.png'), dpi=150)
    plt.close()
    
    # 7. Write Performance Comparison Report
    report_path = os.path.join(args.output_dir, 'ml_performance_comparison_report.md')
    logging.info(f"Writing performance comparison report to {report_path}...")
    with open(report_path, 'w') as fh:
        fh.write("# Comparative Machine Learning CNV Classifier Report\n\n")
        fh.write("This report presents the cross-validation performance metrics of five machine learning classifiers compared to standard concordance-based filtering, evaluated using 95% Confidence Intervals derived from bootstrap resampling.\n\n")
        
        fh.write("## 1. Metrics Comparison Table\n\n")
        fh.write("| Model | Sensitivity (%) | Precision (%) | Specificity (%) | F2-Score | MCC |\n")
        fh.write("| :--- | :---: | :---: | :---: | :---: | :---: |\n")
        
        for name, metrics in model_metrics.items():
            fh.write(f"| **{name}** | "
                     f"{metrics['sens']['mean']*100:.2f} ({metrics['sens']['ci_lower']*100:.1f} - {metrics['sens']['ci_upper']*100:.1f}) | "
                     f"{metrics['prec']['mean']*100:.2f} ({metrics['prec']['ci_lower']*100:.1f} - {metrics['prec']['ci_upper']*100:.1f}) | "
                     f"{metrics['spec']['mean']*100:.2f} ({metrics['spec']['ci_lower']*100:.1f} - {metrics['spec']['ci_upper']*100:.1f}) | "
                     f"{metrics['f2']['mean']:.4f} ({metrics['f2']['ci_lower']:.3f} - {metrics['f2']['ci_upper']:.3f}) | "
                     f"{metrics['mcc']['mean']:.4f} ({metrics['mcc']['ci_lower']:.3f} - {metrics['mcc']['ci_upper']:.3f}) |\n")
                     
        fh.write("\n\n## 2. ROC & PR Curves Plot\n")
        fh.write("The comparative validation curves are plotted below:\n\n")
        fh.write("![Comparative Curves](comparative_performance_curves.png)\n")
        
    # 8. Train best model on full data & run SHAP explainer
    # Find model with best MCC score
    best_model_name = max(
        (name for name in classifiers),
        key=lambda name: model_metrics[name]['mcc']['mean']
    )
    logging.info(f"Best model identified based on MCC: {best_model_name}")
    
    best_clf = classifiers[best_model_name]
    if best_model_name == 'Random Forest':
        best_clf.fit(X_imputed, y)
        explainer = shap.TreeExplainer(best_clf)
        shap_values = explainer(X_imputed)
    else:
        best_clf.fit(X, y)
        explainer = shap.TreeExplainer(best_clf)
        shap_values = explainer(X)
        
    # Plot SHAP summary beeswarm plot
    logging.info("Generating SHAP feature importance beeswarm plot...")
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_imputed if best_model_name == 'Random Forest' else X, show=False)
    plt.title(f"Global Feature Importance (SHAP) - {best_model_name}")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'shap_summary_beeswarm.png'), dpi=150)
    plt.close()
    
    logging.info("ML comparison execution completed successfully!")

if __name__ == '__main__':
    main()
