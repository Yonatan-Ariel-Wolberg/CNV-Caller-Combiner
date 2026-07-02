#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: evaluate_caller_performance.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
"""
Evaluate CNV caller performance (precision and sensitivity) at the probe level.

Truth set and call set BED files must contain 5 tab-separated columns in this
order:
    CHR  START  STOP  CNV_TYPE  SAMPLE_ID

The probes (capture target) BED file must contain at least 3 tab-separated
columns:
    CHR  START  STOP  [...]

Usage:
    evaluate_caller_performance.py \\
        --truth_bed   truthset.bed \\
        --callset_bed callset.bed \\
        --probes_bed  probes.bed \\
        --output      performance.txt
"""

import argparse
import math
import sys

import pandas as pd


def load_bed_file(bed_file):
    """Load a capture-target BED file (CHR, START, STOP[, ...])."""
    return pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end'],
                       usecols=[0, 1, 2])


def load_cnv_file(cnv_file):
    """Load a 5-column CNV BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)."""
    return pd.read_csv(cnv_file, sep='\t', header=None,
                       names=['chr', 'start', 'end', 'cnv_type', 'sample'])


def find_overlaps(probes, cnvs):
    overlaps = []
    for index, probe in probes.iterrows():
        for _, cnv in cnvs.iterrows():
            if (probe['chr'] == cnv['chr'] and
                    probe['start'] < cnv['end'] and
                    probe['end'] > cnv['start']):
                overlaps.append((probe, cnv))
    return overlaps


def categorize_probes(probes, truth_cnv, callset_cnv):
    categorized = {
        'TP': [],
        'FN': [],
        'FP': [],
        'TN': [],
    }

    truth_overlaps = find_overlaps(probes, truth_cnv)
    for probe in probes.itertuples(index=False):
        is_in_truth = any(
            (probe.chr == overlap[0]['chr'] and
             probe.start < overlap[0]['end'] and
             probe.end > overlap[0]['start'])
            for overlap in truth_overlaps
        )

        is_called = any(
            (probe.chr == cnv.chr and
             probe.start < cnv.end and
             probe.end > cnv.start)
            for cnv in callset_cnv.itertuples(index=False)
        )

        if is_in_truth and is_called:
            categorized['TP'].append(probe)
        elif is_in_truth and not is_called:
            categorized['FN'].append(probe)
        elif not is_in_truth and is_called:
            categorized['FP'].append(probe)
        else:
            categorized['TN'].append(probe)

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
