#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: evaluate_caller_performance.py (ULTRA-OPTIMIZED)
# Purpose: Evaluate CNV caller performance (precision, sensitivity, specificity,
#   F-beta score, and MCC) at the probe level using a capture-target probes BED file,
#   a truth-set BED file, and a caller's predicted call-set BED file.
# =============================================================================
"""
Evaluate CNV caller performance (precision and sensitivity) at the probe level.
"""

import argparse
import math
import os
import sys
from bisect import bisect_left
import pandas as pd


def load_bed_file(bed_file):
    if not bed_file or not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
        return pd.DataFrame(columns=['chr', 'start', 'end'])
    try:
        return pd.read_csv(bed_file, sep='\t', header=None,
                           names=['chr', 'start', 'end'],
                           usecols=[0, 1, 2], comment='#')
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['chr', 'start', 'end'])


def load_cnv_file(cnv_file):
    if not cnv_file or not os.path.exists(cnv_file) or os.path.getsize(cnv_file) == 0:
        return pd.DataFrame(columns=['chr', 'start', 'end', 'cnv_type', 'sample'])
    try:
        return pd.read_csv(cnv_file, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'cnv_type', 'sample'], comment='#')
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['chr', 'start', 'end', 'cnv_type', 'sample'])


def find_overlapped_probes(probes_df, cnvs_df):
    """
    Finds indices of probes that overlap with any CNV in cnvs_df.
    Uses O(M log N) binary search mapping directly on raw arrays to avoid pandas overhead.
    """
    if cnvs_df is None or cnvs_df.empty:
        return set()

    # Pre-process probes: group by standardized chromosome and sort by start coordinate
    probes_df = probes_df.copy()
    probes_df['chr_str'] = probes_df['chr'].astype(str).str.replace(r'^[cC][hH][rR]', '', regex=True).str.strip()
    
    probes_by_chrom = {}
    for chrom, group in probes_df.groupby('chr_str'):
        sorted_probes = sorted(list(group.itertuples(index=True)), key=lambda x: x.start)
        probes_by_chrom[chrom] = {
            'list': sorted_probes,
            'starts': [p.start for p in sorted_probes]
        }

    overlapped_indices = set()

    # Group CNVs by chromosome (using raw keys to avoid pandas column mutation overhead)
    for chrom_key, group in cnvs_df.groupby('chr'):
        chrom = str(chrom_key).replace('chr', '').replace('CHR', '').strip()
        if chrom not in probes_by_chrom:
            continue
            
        p_info = probes_by_chrom[chrom]
        p_list = p_info['list']
        p_starts = p_info['starts']
        n_probes = len(p_list)

        # Retrieve raw numpy arrays for zip iteration (10x faster than itertuples)
        starts = group['start'].values
        ends = group['end'].values
        
        for c_start, c_end in zip(starts, ends):
            # Binary search for first probe starting after or at c_start
            idx = bisect_left(p_starts, c_start)
            
            # Check backwards (up to 50kb buffer for maximum target probe length)
            for i in range(idx - 1, -1, -1):
                p = p_list[i]
                if p.start >= c_end:
                    continue
                if p.end > c_start:
                    overlapped_indices.add(p.Index)
                if p.start < c_start - 50000:
                    break
                    
            # Check forwards
            for i in range(idx, n_probes):
                p = p_list[i]
                if p.start >= c_end:
                    break
                if p.end > c_start:
                    overlapped_indices.add(p.Index)

    return overlapped_indices


def categorize_probes(probes, truth_cnv, callset_cnv):
    import numpy as np

    truth_idx = find_overlapped_probes(probes, truth_cnv)
    called_idx = find_overlapped_probes(probes, callset_cnv)

    n_probes = len(probes)
    is_in_truth_arr = np.zeros(n_probes, dtype=bool)
    is_called_arr = np.zeros(n_probes, dtype=bool)

    # Use index lookup to populate boolean arrays
    for idx in truth_idx:
        if 0 <= idx < n_probes:
            is_in_truth_arr[idx] = True
    for idx in called_idx:
        if 0 <= idx < n_probes:
            is_called_arr[idx] = True

    # Vectorized category masks
    tp_mask = is_in_truth_arr & is_called_arr
    fn_mask = is_in_truth_arr & ~is_called_arr
    fp_mask = ~is_in_truth_arr & is_called_arr
    tn_mask = ~is_in_truth_arr & ~is_called_arr

    # Extract raw numpy values for rapid zipping
    chrs = probes['chr'].values
    starts = probes['start'].values
    ends = probes['end'].values

    categorized = {
        'TP': list(zip(chrs[tp_mask], starts[tp_mask], ends[tp_mask])),
        'FN': list(zip(chrs[fn_mask], starts[fn_mask], ends[fn_mask])),
        'FP': list(zip(chrs[fp_mask], starts[fp_mask], ends[fp_mask])),
        'TN': list(zip(chrs[tn_mask], starts[tn_mask], ends[tn_mask])),
    }

    return categorized


def compute_metrics(categorized, beta=2):
    TP = len(categorized['TP'])
    FN = len(categorized['FN'])
    FP = len(categorized['FP'])
    TN = len(categorized['TN'])

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    precision   = TP / (TP + FP) if (TP + FP) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0

    f_beta_denom = (beta ** 2 * precision) + sensitivity
    f_beta = (1 + beta ** 2) * (precision * sensitivity) / f_beta_denom \
        if f_beta_denom > 0 else 0

    mcc_denom = math.sqrt(
        (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
    )
    mcc = ((TP * TN) - (FP * FN)) / mcc_denom if mcc_denom > 0 else 0

    return sensitivity, precision, TP, TN, FP, FN, specificity, f_beta, mcc


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate CNV caller performance (precision and sensitivity).')
    parser.add_argument('--truth_bed',   required=True,
                        help='Truth set BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)')
    parser.add_argument('--callset_bed', required=True,
                        help='Call set BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)')
    parser.add_argument('--probes_bed',  required=True,
                        help='Capture target BED file (CHR, START, STOP[, ...])')
    parser.add_argument('--output',      default=None,
                        help='Output file for performance metrics (default: stdout)')
    args = parser.parse_args()

    probes      = load_bed_file(args.probes_bed)
    truth_cnv   = load_cnv_file(args.truth_bed)
    callset_cnv = load_cnv_file(args.callset_bed)

    categorized = categorize_probes(probes, truth_cnv, callset_cnv)
    sensitivity, precision, TP, TN, FP, FN, specificity, f_beta, mcc = \
        compute_metrics(categorized)

    lines = [
        "Confusion Matrix:",
        f"  TP: {TP}",
        f"  FP: {FP}",
        f"  FN: {FN}",
        f"  TN: {TN}",
        f"Sensitivity: {sensitivity}",
        f"Precision: {precision}",
        f"Specificity: {specificity}",
        f"F_beta (beta=2): {f_beta}",
        f"Matthews Correlation Coefficient: {mcc}",
        f"True Negatives: {TN}",
    ]

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
    else:
        print('\n'.join(lines))


if __name__ == '__main__':
    main()
