# =============================================================================
# Detailed documentation: train_xgboost.py
# Purpose: Trains an XGBoost classifier to distinguish true-positive CNV calls
#   from false-positive noise in combined exome pipeline runs.
# Author/Context: Part of the CNV-Caller-Combiner pipeline.
# Inputs:
#   - --features_dir: Directory containing per-sample feature TSVs
#   - --truth_labels: TSV containing ground truth targets mapping to candidates
#   - --probes_bed: Optional capture-target BED file for probe-overlap matching
# Outputs:
#   - Trained XGBoost model: cnv_model.json
#   - Training metrics report: training_report.txt
#   - ROC and PR curves: roc_curve.svg, pr_curve.svg
#   - SHAP summary plots and values: shap_summary_bar.svg, shap_summary_beeswarm.svg
# Dependencies: XGBoost, Scikit-learn, Imbalanced-learn, Pandas, Numpy, Matplotlib, SHAP
# =============================================================================

import argparse
import glob
import os
import sys


import pandas as pd
import xgboost as xgb
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, precision_recall_curve, auc, average_precision_score
import numpy as np

try:
    from catboost import CatBoostClassifier, Pool
    CATBOOST_AVAILABLE = True
except ImportError:
    CATBOOST_AVAILABLE = False

# All CNV callers supported by the pipeline.  Running all seven callers
# produces fully-populated is_{caller} and qual_norm_{caller} columns in the
# feature matrix.  Callers not provided will have NaN / 0 for their columns.
SUPPORTED_CALLERS = ('canoes', 'clamms', 'xhmm', 'gatk_gcnv', 'cnvkit', 'dragen', 'exomedepth')

# Minimum number of distinct callers required for training the ML classifier.
# With fewer than two callers, the concordance feature is always 1 and
# provides no discriminative signal; at least two callers must contribute to
# the merged input VCF.
MIN_CALLERS_FOR_TRAINING = 2


def _normalise_cnv_type(value):
    """
    Normalizes CNV type labels.

    Args:
        value (str or int): Value representing CNV type (e.g. '1', 'DUP', 'DEL').

    Returns:
        str: 'DUP' for duplications, 'DEL' for deletions, or the raw normalized string.
    """
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    if text in {'1', 'DUP', 'DUPLICATION'}:
        return 'DUP'
    if text in {'0', 'DEL', 'DELETION'}:
        return 'DEL'
    return text


def _overlap_len(a_start, a_end, b_start, b_end):
    """
    Calculates overlap length between two coordinate intervals.

    Args:
        a_start (int): Start of first interval.
        a_end (int): End of first interval.
        b_start (int): Start of second interval.
        b_end (int): End of second interval.

    Returns:
        int: Length of overlapping region (0 if no overlap).
    """
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _load_probes(probes_bed):
    """
    Loads capture-target probes from a BED file.

    Args:
        probes_bed (str): Path to target probes BED file.

    Returns:
        pd.DataFrame: Dataframe with target probes, containing 'probe_id',
                      'chrom', 'start', and 'end'.
    """
    probes = pd.read_csv(
        probes_bed,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end'],
        usecols=[0, 1, 2],
    )
    probes['chrom'] = probes['chrom'].astype(str)
    probes['start'] = probes['start'].astype(int)
    probes['end'] = probes['end'].astype(int)
    probes = probes.reset_index().rename(columns={'index': 'probe_id'})
    return probes


def _probe_ids_for_interval(chrom, start, end, probes_df):
    """
    Finds target probe IDs overlapping a given coordinate interval.

    Args:
        chrom (str): Chromosome name.
        start (int): Start coordinate.
        end (int): End coordinate.
        probes_df (pd.DataFrame): Probes dataframe from _load_probes.

    Returns:
        set: Set of probe IDs that overlap the interval.
    """
    hit = probes_df[
        (probes_df['chrom'] == str(chrom))
        & (probes_df['start'] < int(end))
        & (probes_df['end'] > int(start))
    ]
    return set(hit['probe_id'].tolist())



def merge_features_with_truth_labels(features_df, labels_df, probes_bed=None, min_shared_probes=1):
    """Merge feature rows with truth labels.

    Matching strategy:
      1. Exact key match on sample/chrom/start/end/cnv_type_norm.
      2. Optional fallback for unmatched rows: probe-overlap matching when
         ``probes_bed`` is provided. CNV pairs are considered matches when they
         share at least ``min_shared_probes`` probes on the same sample/chrom/type.
    """
    features = features_df.copy()
    labels = labels_df.copy()
    min_shared_probes = int(min_shared_probes)
    if min_shared_probes < 1:
        raise ValueError('min_shared_probes must be >= 1')

    features['chrom'] = features['chrom'].astype(str)
    labels['chrom'] = labels['chrom'].astype(str)
    features['cnv_type_norm'] = features['cnv_type'].apply(_normalise_cnv_type)
    labels['cnv_type_norm'] = labels['cnv_type'].apply(_normalise_cnv_type)

    labels = labels.reset_index().rename(columns={'index': '_label_row_id'})
    labels_subset = labels[
        ['sample_id', 'chrom', 'start', 'end', 'cnv_type_norm', 'truth_label', '_label_row_id']
    ]
    merged = features.merge(
        labels_subset,
        on=['sample_id', 'chrom', 'start', 'end', 'cnv_type_norm'],
        how='left',
    )

    exact_match_count = int(merged['truth_label'].notna().sum())
    probe_match_count = 0

    if probes_bed:
        unmatched_idx = merged.index[merged['truth_label'].isna()].tolist()
        if unmatched_idx:
            probes = _load_probes(probes_bed)
            used_exact_label_ids = set(
                merged.loc[merged['truth_label'].notna(), '_label_row_id']
                .dropna()
                .astype(int)
                .tolist()
            )
            probe_cache = {}

            # Build candidate matches scored by shared probe count + bp overlap.
            candidates = []
            for f_idx in unmatched_idx:
                feat = merged.loc[f_idx]
                f_key = (feat['chrom'], int(feat['start']), int(feat['end']))
                if f_key not in probe_cache:
                    probe_cache[f_key] = _probe_ids_for_interval(
                        feat['chrom'], feat['start'], feat['end'], probes
                    )
                f_probe_ids = probe_cache[f_key]
                if not f_probe_ids:
                    continue
                lbl_candidates = labels[
                    (labels['sample_id'] == feat['sample_id'])
                    & (labels['chrom'] == feat['chrom'])
                    & (labels['cnv_type_norm'] == feat['cnv_type_norm'])
                    & (~labels['_label_row_id'].isin(used_exact_label_ids))
                ]
                for l_idx, lbl in lbl_candidates.iterrows():
                    l_key = (lbl['chrom'], int(lbl['start']), int(lbl['end']))
                    if l_key not in probe_cache:
                        probe_cache[l_key] = _probe_ids_for_interval(
                            lbl['chrom'], lbl['start'], lbl['end'], probes
                        )
                    l_probe_ids = probe_cache[l_key]
                    shared_probes = len(f_probe_ids & l_probe_ids)
                    if shared_probes < min_shared_probes:
                        continue
                    bp_overlap = _overlap_len(feat['start'], feat['end'], lbl['start'], lbl['end'])
                    candidates.append((
                        shared_probes,
                        bp_overlap,
                        f_idx,
                        int(lbl['_label_row_id']),
                        lbl['truth_label'],
                    ))

            # Greedy one-to-one assignment by best candidate score.
            candidates.sort(reverse=True)
            used_features = set()
            used_labels = set(used_exact_label_ids)
            for _, _, f_idx, label_row_id, truth_label in candidates:
                if f_idx in used_features or label_row_id in used_labels:
                    continue
                merged.at[f_idx, 'truth_label'] = truth_label
                merged.at[f_idx, '_label_row_id'] = label_row_id
                used_features.add(f_idx)
                used_labels.add(label_row_id)
                probe_match_count += 1

    merged = merged[merged['truth_label'].notna()].copy()
    merged['truth_label'] = merged['truth_label'].astype(int)
    merged = merged.drop(columns=['cnv_type_norm', '_label_row_id'], errors='ignore')
    return merged, exact_match_count, probe_match_count


def validate_min_callers(caller_columns):
    """Raise ValueError when fewer than MIN_CALLERS_FOR_TRAINING callers are present.

    Parameters
    ----------
    caller_columns : list[str]
        Names of per-caller flag columns found in the training data.
        These are either ``is_{caller}`` names (feature-extraction output)
        or ``caller_{i}_flag`` names (legacy supp_vec split).

    Raises
    ------
    ValueError
        If ``len(caller_columns) < MIN_CALLERS_FOR_TRAINING``.
    """
    if len(caller_columns) < MIN_CALLERS_FOR_TRAINING:
        raise ValueError(
            f"At least {MIN_CALLERS_FOR_TRAINING} CNV callers are required for "
            f"training the classifier (concordance features are not meaningful "
            f"with fewer callers). "
            f"Found {len(caller_columns)} caller column(s): {list(caller_columns)}. "
            f"Supported callers: {list(SUPPORTED_CALLERS)}."
        )


def prepare_training_data(raw_df, num_callers=7):
    """
    Prepares the features dataframe for ML training.

    Ensures the presence of binary caller indicator features. In legacy workflows,
    parses the SURVIVOR/Truvari support vector (supp_vec) into individual columns.

    Args:
        raw_df (pd.DataFrame): Raw concatenated features.
        num_callers (int): Number of caller columns in the pipeline.

    Returns:
        pd.DataFrame: Formatted features dataframe.
    """
    # Detect named is_* caller columns produced by feature_extraction.py.
    # Exclude non-caller flags such as is_precise and is_sv_refined.
    named_caller_cols = [
        c for c in raw_df.columns
        if c.startswith('is_') and c not in ('is_precise', 'is_sv_refined')
    ]

    if named_caller_cols:
        # Feature-extraction output: is_{caller} columns already present.
        validate_min_callers(named_caller_cols)
        return raw_df

    # Legacy path: raw SURVIVOR/Truvari output with a supp_vec column.
    validate_min_callers([f'caller_{i}_flag' for i in range(num_callers)])
    caller_flags = raw_df['supp_vec'].apply(lambda x: pd.Series(list(str(x))))
    caller_flags.columns = [f'caller_{i}_flag' for i in range(num_callers)]
    caller_flags = caller_flags.astype(int)

    feature_df = pd.concat([raw_df.drop(columns=['supp_vec']), caller_flags], axis=1)
    return feature_df


def train_validation_model(X, y, classifier='catboost'):
    """
    Applies SMOTE oversampling and initializes the classifier (CatBoost or XGBoost).

    Args:
        X (pd.DataFrame): Feature training matrix.
        y (pd.Series or np.ndarray): Ground truth binary labels.
        classifier (str): Classifier architecture ('catboost' or 'xgboost').

    Returns:
        tuple: (model, pd.DataFrame, pd.Series) representing the model
               and resampled features/targets.
    """
    try:
        smote = SMOTE(random_state=42)
        X_res, y_res = smote.fit_resample(X, y)
    except Exception:
        X_res, y_res = X, y

    scale_weight = float(len(y[y == 0])) / float(len(y[y == 1])) if (y == 1).sum() > 0 else 1.0

    if classifier.lower() == 'catboost':
        if not CATBOOST_AVAILABLE:
            raise ImportError("catboost is not installed. Please install catboost or specify --classifier xgboost.")
        model = CatBoostClassifier(
            iterations=500,
            depth=6,
            learning_rate=0.05,
            scale_pos_weight=scale_weight,
            random_state=42,
            verbose=0,
        )
    else:
        model = xgb.XGBClassifier(
            n_estimators=500,
            max_depth=4,
            learning_rate=0.05,
            scale_pos_weight=scale_weight,
            use_label_encoder=False,
            missing=np.nan,
        )

    return model, X_res, y_res


def cross_validate_model(model, X, y, n_splits=5):
    """
    Performs stratified K-Fold cross validation and outputs F1 scores.

    Args:
        model (xgb.XGBClassifier): Classifier to validate.
        X (pd.DataFrame): Input features.
        y (pd.Series or np.ndarray): Ground truth binary labels.
        n_splits (int): Number of validation splits.
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = cross_val_score(model, X, y, cv=skf, scoring='f1')
    
    print(f"Cross-validated F1 scores: {scores}")
    print(f"Mean F1 score: {scores.mean()}")



def _write_line_plot_svg(path, x_values, y_values, title, x_label, y_label, add_diagonal=False):
    width = 800
    height = 600
    margin_left = 80
    margin_right = 40
    margin_top = 60
    margin_bottom = 70
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)
    x_min = float(np.nanmin(x_arr))
    x_max = float(np.nanmax(x_arr))
    y_min = 0.0
    y_max = 1.0

    if x_max == x_min:
        x_max = x_min + 1.0

    def sx(v):
        return margin_left + ((float(v) - x_min) / (x_max - x_min)) * plot_w

    def sy(v):
        return margin_top + (1.0 - ((float(v) - y_min) / (y_max - y_min))) * plot_h

    points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(x_arr, y_arr))
    diag = (
        f'<line x1="{sx(x_min):.2f}" y1="{sy(y_min):.2f}" x2="{sx(x_max):.2f}" '
        f'y2="{sy(y_max):.2f}" stroke="#999" stroke-width="1.5" stroke-dasharray="6 6" />'
        if add_diagonal else ""
    )

    svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white" />
  <text x="{width / 2:.1f}" y="30" text-anchor="middle" font-size="22" font-family="Arial">{title}</text>
  <line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#000" stroke-width="2" />
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height - margin_bottom}" stroke="#000" stroke-width="2" />
  {diag}
  <polyline points="{points}" fill="none" stroke="#1f77b4" stroke-width="3" />
  <text x="{width / 2:.1f}" y="{height - 20}" text-anchor="middle" font-size="16" font-family="Arial">{x_label}</text>
  <text x="20" y="{height / 2:.1f}" text-anchor="middle" font-size="16" font-family="Arial" transform="rotate(-90 20,{height / 2:.1f})">{y_label}</text>
</svg>
"""
    with open(path, 'w') as fh:
        fh.write(svg)


def _write_shap_bar_svg(path, feature_names, mean_abs_shap, top_n=20):
    order = np.argsort(mean_abs_shap)[::-1][:top_n]
    names = [feature_names[i] for i in order]
    vals = [float(mean_abs_shap[i]) for i in order]

    width = 1000
    row_h = 26
    margin_top = 60
    margin_bottom = 30
    margin_left = 320
    margin_right = 60
    plot_w = width - margin_left - margin_right
    height = margin_top + margin_bottom + row_h * max(1, len(names))
    max_v = max(vals) if vals else 1.0

    bars = []
    labels = []
    for idx, (name, value) in enumerate(zip(names, vals)):
        y = margin_top + idx * row_h
        bar_w = (value / max_v) * plot_w if max_v > 0 else 0
        bars.append(
            f'<rect x="{margin_left}" y="{y + 4}" width="{bar_w:.2f}" height="16" fill="#d62728" />'
        )
        labels.append(
            f'<text x="{margin_left - 8}" y="{y + 17}" text-anchor="end" font-size="12" font-family="Arial">{name}</text>'
        )
        labels.append(
            f'<text x="{margin_left + bar_w + 6:.2f}" y="{y + 17}" text-anchor="start" font-size="11" font-family="Arial">{value:.4f}</text>'
        )

    svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white" />
  <text x="{width / 2:.1f}" y="30" text-anchor="middle" font-size="22" font-family="Arial">SHAP summary (mean |value|)</text>
  <line x1="{margin_left}" y1="{margin_top - 10}" x2="{margin_left}" y2="{height - margin_bottom}" stroke="#000" stroke-width="1" />
  {''.join(bars)}
  {''.join(labels)}
</svg>
"""
    with open(path, 'w') as fh:
        fh.write(svg)


def _write_shap_beeswarm_svg(path, feature_names, shap_values, top_n=10, max_points_per_feature=300):
    mean_abs = np.mean(np.abs(shap_values), axis=0)
    order = np.argsort(mean_abs)[::-1][:top_n]
    selected_names = [feature_names[i] for i in order]
    selected_shap = shap_values[:, order]

    width = 1000
    row_h = 40
    margin_top = 70
    margin_bottom = 40
    margin_left = 260
    margin_right = 60
    plot_w = width - margin_left - margin_right
    height = margin_top + margin_bottom + row_h * max(1, len(selected_names))

    x_abs_max = float(np.max(np.abs(selected_shap))) if selected_shap.size else 1.0
    if x_abs_max == 0:
        x_abs_max = 1.0

    def sx(v):
        return margin_left + ((float(v) + x_abs_max) / (2.0 * x_abs_max)) * plot_w

    circles = []
    labels = []
    rng = np.random.default_rng(42)
    for i, name in enumerate(selected_names):
        y_center = margin_top + i * row_h + row_h / 2
        row_vals = selected_shap[:, i]
        if len(row_vals) > max_points_per_feature:
            idx = np.linspace(0, len(row_vals) - 1, max_points_per_feature, dtype=int)
            row_vals = row_vals[idx]
        jitters = rng.uniform(-10, 10, size=len(row_vals))
        for val, jitter in zip(row_vals, jitters):
            circles.append(
                f'<circle cx="{sx(val):.2f}" cy="{(y_center + jitter):.2f}" r="2.5" fill="#1f77b4" fill-opacity="0.55" />'
            )
        labels.append(
            f'<text x="{margin_left - 8}" y="{y_center + 4:.2f}" text-anchor="end" font-size="12" font-family="Arial">{name}</text>'
        )

    x_zero = sx(0.0)
    svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white" />
  <text x="{width / 2:.1f}" y="30" text-anchor="middle" font-size="22" font-family="Arial">SHAP beeswarm (top features)</text>
  <line x1="{x_zero:.2f}" y1="{margin_top - 12}" x2="{x_zero:.2f}" y2="{height - margin_bottom}" stroke="#666" stroke-width="1.5" stroke-dasharray="4 4" />
  {''.join(circles)}
  {''.join(labels)}
  <text x="{width / 2:.1f}" y="{height - 10}" text-anchor="middle" font-size="16" font-family="Arial">SHAP value</text>
</svg>
"""
    with open(path, 'w') as fh:
        fh.write(svg)

def main():
    """Entry point for command-line invocation from Nextflow or the shell.

    Usage
    -----
    python train_xgboost.py \\
        --features_dir  ./out_FEATURES \\
        --truth_labels  truth_labels.tsv \\
        [--output_model  cnv_model.json] \\
        [--output_report training_report.txt]

    The ``features_dir`` is searched recursively for files matching
    ``*_features.tsv`` (the output pattern of ``feature_extraction.py``).

    The ``truth_labels`` TSV must contain at minimum the columns:
        sample_id, chrom, start, end, cnv_type, truth_label
    where ``truth_label`` is 1 for a true CNV and 0 for a false positive.
    """
    parser = argparse.ArgumentParser(
        description='Train an XGBoost classifier on CNV feature matrices.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--features_dir', required=True,
        help='Directory (searched recursively) containing *_features.tsv files '
             'produced by feature_extraction.py.',
    )
    parser.add_argument(
        '--truth_labels', required=True,
        help='TSV file with columns: sample_id, chrom, start, end, cnv_type, truth_label.',
    )
    parser.add_argument(
        '--probes_bed', default=None,
        help='Optional capture-target BED (CHR, START, END). When provided, '
             'unmatched feature/truth rows are additionally matched by shared probes.',
    )
    parser.add_argument(
        '--min_shared_probes', type=int, default=1,
        help='Minimum number of shared probes required for fallback probe-overlap matching.',
    )
    parser.add_argument(
        '--output_model', default='cnv_model.json',
        help='Output path for the trained XGBoost model (XGBoost JSON format).',
    )
    parser.add_argument(
        '--output_report', default='training_report.txt',
        help='Output path for the training metrics report.',
    )
    parser.add_argument(
        '--output_roc_plot', default='roc_curve.svg',
        help='Output path for ROC curve plot (SVG).',
    )
    parser.add_argument(
        '--output_pr_plot', default='pr_curve.svg',
        help='Output path for precision-recall curve plot (SVG).',
    )
    parser.add_argument(
        '--output_roc_data', default='roc_curve.tsv',
        help='Output path for ROC curve points TSV.',
    )
    parser.add_argument(
        '--output_pr_data', default='pr_curve.tsv',
        help='Output path for precision-recall curve points TSV.',
    )
    parser.add_argument(
        '--output_shap_values', default='shap_values.tsv',
        help='Output path for per-sample SHAP values TSV.',
    )
    parser.add_argument(
        '--output_shap_summary_plot', default='shap_summary_bar.svg',
        help='Output path for SHAP summary bar plot (SVG).',
    )
    parser.add_argument(
        '--output_shap_beeswarm_plot', default='shap_summary_beeswarm.svg',
        help='Output path for SHAP beeswarm plot (SVG).',
    )
    parser.add_argument(
        '--classifier', default='catboost', choices=['catboost', 'xgboost'],
        help='Machine learning classifier architecture to train (catboost or xgboost; default: catboost).',
    )
    parser.add_argument(
        '--shap_top_features', type=int, default=20,
        help='Number of top features to include in SHAP summary plots.',
    )
    args = parser.parse_args()
    if args.min_shared_probes < 1:
        sys.exit('ERROR: --min_shared_probes must be >= 1')

    # ── 1. Load all feature TSV files ─────────────────────────────────────────
    tsv_files = sorted(
        glob.glob(os.path.join(args.features_dir, '**', '*_features.tsv'), recursive=True)
    )
    if not tsv_files:
        tsv_files = sorted(glob.glob(os.path.join(args.features_dir, '*_features.tsv')))

    if not tsv_files:
        sys.exit(
            f"ERROR: No '*_features.tsv' files found in '{args.features_dir}'. "
            "Ensure the feature_extraction workflow completed successfully."
        )

    dfs = [pd.read_csv(f, sep='\t') for f in tsv_files]
    features_df = pd.concat(dfs, ignore_index=True)

    # ── 2. Load and validate truth labels ────────────────────────────────────
    labels_df = pd.read_csv(args.truth_labels, sep='\t')
    required_cols = {'sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label'}
    missing_cols = required_cols - set(labels_df.columns)
    if missing_cols:
        sys.exit(
            f"ERROR: truth_labels file is missing required columns: {missing_cols}"
        )

    # ── 3. Merge features with truth labels ──────────────────────────────────
    merged, exact_match_count, probe_match_count = merge_features_with_truth_labels(
        features_df,
        labels_df[['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label']],
        probes_bed=args.probes_bed,
        min_shared_probes=args.min_shared_probes,
    )

    if merged.empty:
        sys.exit(
            "ERROR: No matching records found after joining features with truth labels. "
            "Check that sample_id / chrom / start / end / cnv_type values align. "
            "If CNV coordinates differ, provide --probes_bed to enable probe-overlap matching."
        )

    y = merged['truth_label']
    X_raw = merged.drop(columns=['truth_label'])

    # prepare_training_data validates caller columns and handles legacy supp_vec
    X_raw = prepare_training_data(X_raw)

    # Drop non-numeric, identifier, and capture-biased breakpoint flags before passing to the model.
    # Note: 'is_precise' and 'is_sv_refined' are intentionally excluded because breakpoint refinement
    # via split-reads/discordant pairs is inherently limited in short-read exome sequencing, where
    # breakpoints fall in uncaptured intronic space. Penalizing unrefined calls would cause false-negative
    # dropouts of genuine macro-deletions and multi-exon CNVs.
    excluded_cols = [c for c in ('sample_id', 'chrom', 'cnv_type', 'is_precise', 'is_sv_refined') if c in X_raw.columns]
    X = X_raw.drop(columns=excluded_cols).select_dtypes(include=[np.number])

    # ── 4. Train and cross-validate model ─────────────────────────────────────
    model, X_res, y_res = train_validation_model(X, y, classifier=args.classifier)
    cross_validate_model(model, X_res, y_res)
    model.fit(X_res, y_res)

    # ── 5. Persist the trained model ──────────────────────────────────────────
    if hasattr(model, 'save_model'):
        try:
            model.save_model(args.output_model, format="json" if args.output_model.endswith(".json") else "cbm")
        except TypeError:
            model.save_model(args.output_model)

    # ── 6. Evaluation curves (ROC / PR) ──────────────────────────────────────
    y_score = model.predict_proba(X)[:, 1]
    fpr, tpr, roc_thresholds = roc_curve(y, y_score)
    precision, recall, pr_thresholds = precision_recall_curve(y, y_score)
    roc_auc = auc(fpr, tpr)
    pr_auc = average_precision_score(y, y_score)

    pd.DataFrame({
        'fpr': fpr,
        'tpr': tpr,
        'threshold': roc_thresholds,
    }).to_csv(args.output_roc_data, sep='\t', index=False)

    pr_thresholds_padded = np.append(pr_thresholds, np.nan)
    pd.DataFrame({
        'recall': recall,
        'precision': precision,
        'threshold': pr_thresholds_padded,
    }).to_csv(args.output_pr_data, sep='\t', index=False)

    _write_line_plot_svg(
        args.output_roc_plot,
        fpr,
        tpr,
        title=f'ROC curve (AUC={roc_auc:.4f})',
        x_label='False positive rate',
        y_label='True positive rate',
        add_diagonal=True,
    )
    _write_line_plot_svg(
        args.output_pr_plot,
        recall,
        precision,
        title=f'Precision-Recall curve (AUPRC={pr_auc:.4f})',
        x_label='Recall',
        y_label='Precision',
    )

    # ── 7. SHAP values computation and plotting ──────────────────────────────
    if CATBOOST_AVAILABLE and isinstance(model, CatBoostClassifier):
        pool = Pool(X)
        shap_all = model.get_feature_importance(pool, type='ShapValues')
        shap_values = shap_all[:, :-1]
    else:
        dmatrix = xgb.DMatrix(X, feature_names=list(X.columns))
        contribs = model.get_booster().predict(dmatrix, pred_contribs=True)
        shap_values = contribs[:, :-1]  # last column is the model bias term

    shap_df = pd.DataFrame(shap_values, columns=X.columns)
    shap_df.insert(0, 'sample_index', np.arange(len(shap_df)))
    shap_df.to_csv(args.output_shap_values, sep='\t', index=False)

    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
    _write_shap_bar_svg(
        args.output_shap_summary_plot,
        list(X.columns),
        mean_abs_shap,
        top_n=max(1, args.shap_top_features),
    )
    _write_shap_beeswarm_svg(
        args.output_shap_beeswarm_plot,
        list(X.columns),
        shap_values,
        top_n=max(1, args.shap_top_features),
    )

    # ── 8. Write training report ──────────────────────────────────────────────
    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())
    with open(args.output_report, 'w') as fh:
        fh.write(f"Classifier:                {args.classifier}\n")
        fh.write(f"Feature TSV files:         {len(tsv_files)}\n")
        fh.write(f"Training samples (calls):  {len(X)}\n")
        fh.write(f"Exact key matches:         {exact_match_count}\n")
        fh.write(f"Probe-overlap matches:     {probe_match_count}\n")
        fh.write(f"True CNVs  (label=1):      {n_pos}\n")
        fh.write(f"False CNVs (label=0):      {n_neg}\n")
        fh.write(f"ROC AUC:                   {roc_auc:.6f}\n")
        fh.write(f"PR AUC:                    {pr_auc:.6f}\n")
        fh.write(f"Feature columns:           {list(X.columns)}\n")
        fh.write(f"Model saved to:            {args.output_model}\n")
        fh.write(f"ROC plot:                  {args.output_roc_plot}\n")
        fh.write(f"PR plot:                   {args.output_pr_plot}\n")
        fh.write(f"ROC data:                  {args.output_roc_data}\n")
        fh.write(f"PR data:                   {args.output_pr_data}\n")
        fh.write(f"SHAP values:               {args.output_shap_values}\n")
        fh.write(f"SHAP summary plot:         {args.output_shap_summary_plot}\n")
        fh.write(f"SHAP beeswarm plot:        {args.output_shap_beeswarm_plot}\n")

    print(f"Model saved to {args.output_model}")
    print(f"Report written to {args.output_report}")


if __name__ == '__main__':
    main()

