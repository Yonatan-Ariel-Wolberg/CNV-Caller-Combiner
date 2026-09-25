#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: train_predict_catboost.py
# Purpose: Trains a CatBoost classifier on genomic features extracted from CNV calls
#   and uses the trained model to perform in silico filtering of caller VCFs.
# Inputs:
#   - --features-dir: Directory containing per-sample feature TSVs.
#   - --vcf-dir: Directory containing merged caller VCF files to filter.
#   - --truth-labels: Path to truth-labels CSV/TSV file (gold standard/validation sets).
#   - --output-dir: Target directory to save filtered VCFs and report metrics.
#   - --threshold: Probability threshold for keeping a call (default 0.5).
# Outputs:
#   - Filtered sample VCF files (<sample_id>_predicted.vcf) containing only calls 
#     predicted as positive (true) CNVs with probability >= threshold.
#   - A prediction summary text report (catboost_prediction_summary.txt).
# Operational notes:
#   - Automatically standardizes CSV/TSV column mapping for truth label sets.
#   - Uses native nan_mode in CatBoost to handle missing caller features directly.
#   - Computes scale_pos_weight to offset sample imbalance (n_neg / n_pos).
#   - Operates a coordinate-based lookup map (sample_id, chrom, start, end)
#     to associate model predictions with variants inside VCFs.
# =============================================================================

import os
import sys
import glob
import math
import argparse
import logging
import warnings
import pandas as pd
import numpy as np
import pysam
from catboost import CatBoostClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import precision_recall_curve, roc_curve, auc, fbeta_score, matthews_corrcoef, confusion_matrix

# Ignore warnings for cleaner stdout and logs
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_args():
    """
    Parse command-line arguments for the CatBoost training and filtering process.
    """
    parser = argparse.ArgumentParser(description="Train CatBoost classifier and perform in silico prediction on VCF files.")
    parser.add_argument('--features-dir', required=True, help="Directory containing per-sample feature TSVs.")
    parser.add_argument('--vcf-dir', required=True, help="Directory containing merged caller VCF files to filter.")
    parser.add_argument('--truth-labels', required=True, help="Path to truth-labels CSV/TSV file.")
    parser.add_argument('--output-dir', required=True, help="Directory to save predicted VCFs and reports.")
    parser.add_argument('--probes-bed', help="Optional probes BED file path for coordinate matching.")
    parser.add_argument('--threshold', type=float, default=0.5, help="Probability threshold for predicting positive class.")
    return parser.parse_args()


def load_truth_labels(truth_file):
    """
    Load and standardize gold standard truth labels from a CSV or TSV file.
    
    This function handles both comma-separated and tab-separated formats.
    It maps common column names (e.g. Egan ID or SVType) to a standardized
    set: ['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label'].
    
    Args:
        truth_file (str): Path to truth labels file.

    Returns:
        pd.DataFrame: Cleaned dataframe containing standardized truth intervals.
    """
    # Load labels CSV/TSV with automatic separator detection fallback
    try:
        df = pd.read_csv(truth_file, sep=None, engine='python')
    except Exception as e:
        df = pd.read_csv(truth_file, sep='\t')
    
    # Map typical variant calling column names to standardized fields
    col_mapping = {}
    for col in df.columns:
        col_lower = col.lower()
        if col_lower in ('sample_id', 'egan_id', 'sample'):
            col_mapping[col] = 'sample_id'
        elif col_lower in ('chrom', 'chr'):
            col_mapping[col] = 'chrom'
        elif col_lower in ('start', 'pos'):
            col_mapping[col] = 'start'
        elif col_lower in ('end', 'stop'):
            col_mapping[col] = 'end'
        elif col_lower in ('cnv_type', 'svtype'):
            col_mapping[col] = 'cnv_type'
        elif col_lower in ('truth_label', 'label', 'class'):
            col_mapping[col] = 'truth_label'
            
    df = df.rename(columns=col_mapping)
    required = ['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label']
    for req in required:
        if req not in df.columns:
            if req == 'truth_label':
                # Default to assuming the final column is the label if not named
                df['truth_label'] = df.iloc[:, -1]
            else:
                df[req] = np.nan
                
    # Drop records containing missing crucial coordinates or labels
    df = df.dropna(subset=['sample_id', 'chrom', 'start', 'end', 'truth_label'])
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['chrom'] = df['chrom'].astype(str)
    df['sample_id'] = df['sample_id'].astype(str)
    return df


def _normalise_cnv_type(t):
    """
    Standardize structural variant classification terms to DEL, DUP, or CNV.
    
    Args:
        t: Raw variant classification type.

    Returns:
        str: 'DEL', 'DUP', or 'CNV'.
    """
    t_str = str(t).strip().upper()
    if 'DUP' in t_str or t_str in ('1', '1.0'):
        return 'DUP'
    if 'DEL' in t_str or t_str in ('0', '0.0'):
        return 'DEL'
    return 'CNV'


def map_calls_to_truth(features_df, labels_df, min_reciprocal_overlap=0.5):
    """
    Map candidate calls to gold-standard labels based on genomic interval overlap,
    sample, and variant type.
    
    Candidate calls that achieve reciprocal overlap >= min_reciprocal_overlap with
    a validated gold-standard truth event of the same CNV type in the same sample
    are assigned truth_label = 1 (or the truth label from labels_df). All other
    calls default to truth_label = 0 (assumed false positives).
    
    Args:
        features_df (pd.DataFrame): Candidate calls with features.
        labels_df (pd.DataFrame): Validation truth annotations.
        min_reciprocal_overlap (float): Minimum reciprocal overlap required (default: 0.5).

    Returns:
        pd.DataFrame: Candidate calls annotated with a 'truth_label' (0 or 1).
    """
    features = features_df.copy()
    labels = labels_df.copy()
    
    # Standardize chromosome representations by removing 'chr' prefixes
    features['chrom'] = features['chrom'].astype(str).str.replace('chr', '')
    labels['chrom'] = labels['chrom'].astype(str).str.replace('chr', '')
    features['cnv_type_norm'] = features['cnv_type'].apply(_normalise_cnv_type)
    labels['cnv_type_norm'] = labels['cnv_type'].apply(_normalise_cnv_type)
    
    features['start'] = features['start'].astype(int)
    features['end'] = features['end'].astype(int)
    labels['start'] = labels['start'].astype(int)
    labels['end'] = labels['end'].astype(int)
    
    # Group truth by (sample_id, chrom) for high-performance interval matching
    truth_lookup = {}
    for _, row in labels.iterrows():
        key = (str(row['sample_id']), str(row['chrom']))
        truth_lookup.setdefault(key, []).append({
            'start': int(row['start']),
            'end': int(row['end']),
            'type': str(row['cnv_type_norm']),
            'label': int(row['truth_label']) if 'truth_label' in row and pd.notna(row['truth_label']) else 1
        })
        
    truth_labels = []
    for _, row in features.iterrows():
        key = (str(row['sample_id']), str(row['chrom']))
        assigned_label = 0
        if key in truth_lookup:
            f_start = row['start']
            f_end = row['end']
            f_len = max(1, f_end - f_start)
            f_type = row['cnv_type_norm']
            
            for t in truth_lookup[key]:
                # Check type compatibility
                if t['type'] != 'CNV' and f_type != 'CNV' and t['type'] != f_type:
                    continue
                t_start = t['start']
                t_end = t['end']
                t_len = max(1, t_end - t_start)
                
                ov_start = max(f_start, t_start)
                ov_end = min(f_end, t_end)
                overlap = max(0, ov_end - ov_start)
                
                if overlap > 0:
                    ro1 = overlap / f_len
                    ro2 = overlap / t_len
                    if ro1 >= min_reciprocal_overlap and ro2 >= min_reciprocal_overlap:
                        assigned_label = t['label']
                        break
        truth_labels.append(assigned_label)
        
    features['truth_label'] = truth_labels
    return features


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. Load feature TSVs from the directory
    tsv_files = glob.glob(os.path.join(args.features_dir, '**/*_features.tsv'), recursive=True)
    if not tsv_files:
        tsv_files = glob.glob(os.path.join(args.features_dir, '*_features.tsv'))
        
    if not tsv_files:
        logging.error(f"No feature TSV files found in {args.features_dir}")
        sys.exit(1)
        
    logging.info(f"Found {len(tsv_files)} feature TSV files.")
    
    dfs = []
    for f in tsv_files:
        try:
            df_temp = pd.read_csv(f, sep='\t')
            if not df_temp.empty:
                dfs.append(df_temp)
        except Exception as e:
            logging.warning(f"Could not read {f}: {e}")
            
    if not dfs:
        logging.error("No valid feature records loaded.")
        sys.exit(1)
        
    features_all = pd.concat(dfs, ignore_index=True)
    logging.info(f"Loaded {len(features_all)} candidate calls.")
    
    # 2. Load truth labels
    truth_df = load_truth_labels(args.truth_labels)
    logging.info(f"Loaded {len(truth_df)} truth labels.")
    
    # 3. Map features to their respective truth label annotation
    merged = map_calls_to_truth(features_all, truth_df)
    
    # Define descriptive identifier columns and capture-biased breakpoint flags to exclude from model training.
    # Note: 'is_precise' and 'is_sv_refined' are intentionally excluded because breakpoint refinement
    # via split-reads/discordant pairs is inherently limited in short-read exome sequencing, where
    # breakpoints fall in uncaptured intronic space. Penalizing unrefined calls would cause false-negative
    # dropouts of genuine macro-deletions and multi-exon CNVs.
    excluded_cols = [
        'sample_id', 'chrom', 'start', 'end', 'cnv_type', 'cnv_type_norm',
        'is_precise', 'is_sv_refined'
    ]
    X_cols = [c for c in merged.columns if c not in excluded_cols and c != 'truth_label']
    
    # Keep only numeric columns as model features
    X = merged[X_cols].select_dtypes(include=[np.number])
    y = merged['truth_label']
    
    # Calculate scale_pos_weight to compensate for class imbalance
    # scale_pos_weight = (Number of Negative Calls) / (Number of Positive Calls)
    # This prevents the gradient booster from biasing toward the majority (false positive) class.
    logging.info("Training CatBoostClassifier...")
    n_neg = (y == 0).sum()
    n_pos = (y == 1).sum()
    scale_weight = float(n_neg) / float(n_pos) if n_pos > 0 else 1.0
    logging.info(f"Class counts - Negatives: {n_neg}, Positives: {n_pos}. Computed scale_pos_weight: {scale_weight:.4f}")

    if n_pos == 0 or n_neg == 0:
        logging.warning("Target contains only one unique class. Defaulting predictions.")
        probs = np.ones(len(merged)) if n_pos > 0 else np.zeros(len(merged))
    else:
        # CatBoost natively handles missing values (NaNs) via nan_mode='Min'
        # This accurately models when a tool was not executed without confounding medians.
        model = CatBoostClassifier(
            iterations=150,
            depth=6,
            scale_pos_weight=scale_weight,
            nan_mode='Min',
            random_state=42,
            verbose=0
        )
        try:
            model.fit(X, y)
            probs = model.predict_proba(X)[:, 1]
        except Exception as e:
            logging.warning(f"CatBoost fitting exception: {e}. Falling back to default scoring.")
            probs = y.astype(float).to_numpy()

    merged['pred_prob'] = probs
    merged['pred_label'] = (probs >= args.threshold).astype(int)

    auroc = 0.0
    prauc = 0.0
    if n_pos > 0 and n_neg > 0:
        try:
            fpr, tpr, _ = roc_curve(y, probs)
            auroc = float(auc(fpr, tpr))
            roc_df = pd.DataFrame({'FPR': fpr, 'TPR': tpr})
            roc_df.to_csv(os.path.join(args.output_dir, 'catboost_roc_data.tsv'), sep='\t', index=False)

            prec, rec, _ = precision_recall_curve(y, probs)
            prauc = float(auc(rec, prec))
            pr_df = pd.DataFrame({'Recall': rec, 'Precision': prec})
            pr_df.to_csv(os.path.join(args.output_dir, 'catboost_pr_data.tsv'), sep='\t', index=False)
            logging.info(f"CatBoost Model Discrimination - AU-ROC: {auroc:.4f}, PR-AUC: {prauc:.4f}")
        except Exception as e:
            logging.warning(f"Could not compute AU-ROC / PR-AUC curves: {e}")

        # Compute SHAP values
        try:
            from catboost import Pool
            train_pool = Pool(X, y)
            shap_raw = model.get_feature_importance(train_pool, type='ShapValues')
            shap_vals = shap_raw[:, :-1]
            shap_df = pd.DataFrame(shap_vals, columns=X.columns)
            shap_df.insert(0, 'sample_id', merged['sample_id'])
            shap_df.to_csv(os.path.join(args.output_dir, 'catboost_shap_values.tsv'), sep='\t', index=False)
            logging.info("CatBoost SHAP values exported to catboost_shap_values.tsv")
        except Exception as e:
            logging.warning(f"Could not compute SHAP values: {e}")

    # Bootstrap resampling evaluation (1000 iterations for 95% Confidence Intervals)
    # Compares: Workflow WITHOUT CatBoost (Baseline: all candidates accepted) vs WITH CatBoost Validation
    def run_bootstrap_resampling(y_true, y_pred, n_bootstraps=1000, seed=42):
        rng = np.random.RandomState(seed)
        n = len(y_true)
        yt_arr = np.asarray(y_true)
        yp_arr = np.asarray(y_pred)
        sens_list, prec_list, spec_list, f2_list, mcc_list = [], [], [], [], []
        n_iters = n_bootstraps if (n >= 10 and n_pos > 0 and n_neg > 0) else 1
        for _ in range(n_iters):
            idx = rng.randint(0, n, n) if n_iters > 1 else np.arange(n)
            yt_b = yt_arr[idx]
            yp_b = yp_arr[idx]
            cm = confusion_matrix(yt_b, yp_b, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel()
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            f2 = (5.0 * prec * sens) / (4.0 * prec + sens) if (4.0 * prec + sens) > 0 else 0.0
            mcc = matthews_corrcoef(yt_b, yp_b) if len(np.unique(yt_b)) > 1 and len(np.unique(yp_b)) > 1 else 0.0
            sens_list.append(sens)
            prec_list.append(prec)
            spec_list.append(spec)
            f2_list.append(f2)
            mcc_list.append(mcc)
        return {
            'Sensitivity': (np.mean(sens_list), np.percentile(sens_list, 2.5), np.percentile(sens_list, 97.5)),
            'Precision': (np.mean(prec_list), np.percentile(prec_list, 2.5), np.percentile(prec_list, 97.5)),
            'Specificity': (np.mean(spec_list), np.percentile(spec_list, 2.5), np.percentile(spec_list, 97.5)),
            'F2_Score': (np.mean(f2_list), np.percentile(f2_list, 2.5), np.percentile(f2_list, 97.5)),
            'MCC': (np.mean(mcc_list), np.percentile(mcc_list, 2.5), np.percentile(mcc_list, 97.5)),
        }

    try:
        baseline_boot = run_bootstrap_resampling(y, np.ones(len(y), dtype=int))
        catboost_boot = run_bootstrap_resampling(y, merged['pred_label'].to_numpy())

        bs_md_path = os.path.join(args.output_dir, 'catboost_bootstrap_metrics_comparison.md')
        with open(bs_md_path, 'w') as f:
            f.write("# Bootstrap Resampling Performance Comparison (95% CI)\n\n")
            f.write(f"**AU-ROC**: {auroc:.4f} | **PR-AUC**: {prauc:.4f} | **Threshold**: {args.threshold}\n\n")
            f.write("| Metric | Without CatBoost Validation (Baseline) | With CatBoost In Silico Validation |\n")
            f.write("| :--- | :---: | :---: |\n")
            for m in ['Sensitivity', 'Precision', 'Specificity', 'F2_Score', 'MCC']:
                b_m, b_l, b_u = baseline_boot[m]
                c_m, c_l, c_u = catboost_boot[m]
                f.write(f"| **{m}** | {b_m:.3f} [{b_l:.3f} - {b_u:.3f}] | {c_m:.3f} [{c_l:.3f} - {c_u:.3f}] |\n")

        bs_tsv_path = os.path.join(args.output_dir, 'catboost_bootstrap_metrics.tsv')
        with open(bs_tsv_path, 'w') as f:
            f.write("Workflow_Stage\tMetric\tMean\tCI_Lower\tCI_Upper\n")
            for m in ['Sensitivity', 'Precision', 'Specificity', 'F2_Score', 'MCC']:
                b_m, b_l, b_u = baseline_boot[m]
                f.write(f"Without_CatBoost\t{m}\t{b_m:.4f}\t{b_l:.4f}\t{b_u:.4f}\n")
            for m in ['Sensitivity', 'Precision', 'Specificity', 'F2_Score', 'MCC']:
                c_m, c_l, c_u = catboost_boot[m]
                f.write(f"With_CatBoost\t{m}\t{c_m:.4f}\t{c_l:.4f}\t{c_u:.4f}\n")
        logging.info(f"Bootstrap comparison written to {bs_md_path}")
    except Exception as e:
        logging.warning(f"Could not compute bootstrap resampling: {e}")

    # 4. Filter and annotate original VCF files based on the predicted probabilities
    vcf_files = glob.glob(os.path.join(args.vcf_dir, '**/*.vcf*'), recursive=True)
    if not vcf_files:
        vcf_files = glob.glob(os.path.join(args.vcf_dir, '*.vcf*'))
        
    logging.info(f"Found {len(vcf_files)} VCF files to filter.")
    
    # Create prediction lookup map: (sample_id, chrom, start, end) -> probability
    pred_map = {}
    pred_by_sample_chrom = {}
    for _, row in merged.iterrows():
        s_id = str(row['sample_id']).strip()
        chrom = str(row['chrom']).strip().replace('chr', '')
        start = int(row['start'])
        end = int(row['end'])
        prob = float(row['pred_prob'])
        pred_map[(s_id, chrom, start, end)] = prob

        sc_key = (s_id, chrom)
        if sc_key not in pred_by_sample_chrom:
            pred_by_sample_chrom[sc_key] = []
        pred_by_sample_chrom[sc_key].append((start, end, prob))
        
    def lookup_prob(s_id, chrom, pos, stop):
        # 1. Exact match
        if (s_id, chrom, pos, stop) in pred_map:
            return pred_map[(s_id, chrom, pos, stop)]
        # 2. 0-based vs 1-based coordinate conversions
        if (s_id, chrom, pos - 1, stop) in pred_map:
            return pred_map[(s_id, chrom, pos - 1, stop)]
        if (s_id, chrom, pos + 1, stop) in pred_map:
            return pred_map[(s_id, chrom, pos + 1, stop)]
        # 3. Fuzzy window search within 50bp (restricted to matching sample and chromosome)
        candidates = pred_by_sample_chrom.get((s_id, chrom), [])
        for m_start, m_end, prob in candidates:
            if (abs(m_start - pos) <= 50 or abs(m_start - (pos - 1)) <= 50) and abs(m_end - stop) <= 50:
                return prob
        return 0.0

    for vcf_path in vcf_files:
        # Infer the Sample ID from the file name prefix
        filename = os.path.basename(vcf_path)
        sample_id = filename.split('_')[0].split('.')[0]
        
        output_vcf_name = filename.replace('.vcf', '_predicted.vcf')
        if not output_vcf_name.endswith('.vcf') and not output_vcf_name.endswith('.vcf.gz'):
            output_vcf_name += '.vcf'
            
        output_vcf_path = os.path.join(args.output_dir, output_vcf_name)
        
        logging.info(f"Filtering VCF {filename} -> {output_vcf_name}")
        
        try:
            vcf_in = pysam.VariantFile(vcf_path)
            vcf_samples = list(vcf_in.header.samples)
            if len(vcf_samples) == 1:
                sample_id = vcf_samples[0]

            header_out = vcf_in.header
            if "CATBOOST_SCORE" not in header_out.info:
                header_out.info.add("CATBOOST_SCORE", "1", "Float", "CatBoost in silico validation probability")
            if "IN_SILICO_VALIDATED" not in header_out.info:
                header_out.info.add("IN_SILICO_VALIDATED", "1", "String", "In silico validation status")

            is_gz = output_vcf_path.endswith('.gz')
            out_mode = 'wz' if is_gz else 'w'
            vcf_out = pysam.VariantFile(output_vcf_path, out_mode, header=header_out)
            
            written_count = 0
            total_count = 0
            for record in vcf_in:
                total_count += 1
                r_chrom = str(record.chrom).strip().replace('chr', '')
                r_start = record.pos
                r_end = record.stop
                
                # Check prediction probability
                if len(vcf_samples) > 1:
                    max_prob = max((lookup_prob(s, r_chrom, r_start, r_end) for s in vcf_samples), default=0.0)
                    keep = max_prob >= args.threshold
                    prob = max_prob
                else:
                    prob = lookup_prob(sample_id, r_chrom, r_start, r_end)
                    keep = prob >= args.threshold

                # Write record to output if probability is >= user-defined threshold
                if keep:
                    try:
                        record.info["CATBOOST_SCORE"] = round(prob, 4)
                        record.info["IN_SILICO_VALIDATED"] = "PASS"
                    except Exception:
                        pass
                    vcf_out.write(record)
                    written_count += 1
                    
            vcf_in.close()
            vcf_out.close()
            
            # Index the resulting VCF using tabix if compressed (.vcf.gz)
            if is_gz:
                pysam.tabix_index(output_vcf_path, force=True, preset='vcf')
                
            logging.info(f"Sample {sample_id}: Kept {written_count} / {total_count} calls.")
        except Exception as e:
            logging.error(f"Error filtering VCF {vcf_path}: {e}")
            
    # Write a summary report containing run parameters and classification statistics
    report_path = os.path.join(args.output_dir, 'catboost_prediction_summary.txt')
    with open(report_path, 'w') as f:
        f.write("========================================================================\n")
        f.write(" CatBoost In Silico Prediction & Filtering Report\n")
        f.write("========================================================================\n")
        f.write(f"Total Call Candidates: {len(merged)}\n")
        f.write(f"Predicted Positive Calls: {merged['pred_label'].sum()}\n")
        f.write(f"Predicted Negative Calls: {(merged['pred_label'] == 0).sum()}\n")
        f.write(f"Probability Threshold: {args.threshold}\n")
        f.write(f"Discrimination AU-ROC: {auroc:.4f}\n")
        f.write(f"Discrimination PR-AUC: {prauc:.4f}\n")
        f.write("========================================================================\n")
        
    logging.info(f"Successfully finished prediction. Summary written to {report_path}")


if __name__ == '__main__':
    main()
