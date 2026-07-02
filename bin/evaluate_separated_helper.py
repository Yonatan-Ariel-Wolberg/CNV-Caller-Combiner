#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: evaluate_separated_helper.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import argparse
import bisect
import math
import sys
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate CNV caller performance on AnnotSV ACMG score-separated CSVs.")
    parser.add_argument('--input-csv', required=True, help="Path to the score-separated CNVs CSV file.")
    parser.add_argument('--map-file', required=True, help="Path to the Sample ID .map file (EGAF to EGAN).")
    parser.add_argument('--truth-bed', required=True, help="Path to the truthset BED file.")
    parser.add_argument('--probes-bed', required=True, help="Path to the target probes BED file.")
    parser.add_argument('--workflow-type', required=True, help="Workflow type (e.g. Control or Mosdepth).")
    parser.add_argument('--output', required=True, help="Path to write the performance report.")
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

def load_probes_bed(probes_file):
    """Load target probes BED file (CHR, START, STOP)."""
    df = pd.read_csv(probes_file, sep='\t', header=None, comment='#')
    # Standardize to 3 columns
    df = df.iloc[:, [0, 1, 2]]
    df.columns = ['chr', 'start', 'end']
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['chr'] = df['chr'].astype(str)
    return df

def load_truth_bed(truth_file):
    """Load 5-column truth BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)."""
    df = pd.read_csv(truth_file, sep='\t', header=None, comment='#', keep_default_na=False)
    # Standardize to 5 columns
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
    # Filter out empty samples
    df = df[df['sample'] != ""]
    return df

def find_overlapping_probe_indices(probes, cnvs):
    """Optimized sweep algorithm to find indices of probes overlapping with CNVs."""
    if cnvs.empty:
        return set()
        
    overlapping_indices = set()
    
    # Sort and group by chromosome to make search O(N log M)
    probes_by_chr = {c: group for c, group in probes.groupby('chr')}
    cnvs_by_chr = {c: group.sort_values('start') for c, group in cnvs.groupby('chr')}
    
    for chrom, p_group in probes_by_chr.items():
        if chrom not in cnvs_by_chr:
            continue
        c_group = cnvs_by_chr[chrom]
        c_starts = c_group['start'].values
        c_ends = c_group['end'].values
        
        p_starts = p_group['start'].values
        p_ends = p_group['end'].values
        p_indices = p_group.index.values
        
        for idx_in_group in range(len(p_group)):
            p_start = p_starts[idx_in_group]
            p_end = p_ends[idx_in_group]
            p_idx = p_indices[idx_in_group]
            
            # Find CNVs that start before p_end
            limit = bisect.bisect_left(c_starts, p_end)
            for i in range(limit):
                if c_ends[i] > p_start:
                    overlapping_indices.add(p_idx)
                    break # A probe only needs to overlap once
                    
    return overlapping_indices

def compute_metrics(tp, fn, fp, tn, beta=2):
    """Compute sensitivity, precision, specificity, F_beta, and MCC."""
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    f_beta_denom = (beta ** 2 * precision) + sensitivity
    f_beta = (1 + beta ** 2) * (precision * sensitivity) / f_beta_denom if f_beta_denom > 0 else 0.0

    mcc_denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / mcc_denom if mcc_denom > 0 else 0.0

    return sensitivity, precision, specificity, f_beta, mcc

def main():
    args = parse_args()
    
    # 1. Load Sample ID mapping
    egaf_to_egan = load_map_file(args.map_file)
    cohort_egan_samples = list(egaf_to_egan.values())
    
    # 2. Load callset CNVs
    try:
        callset_csv = pd.read_csv(args.input_csv, keep_default_na=False)
    except Exception as e:
        logging.error(f"Failed to read input CSV: {e}")
        sys.exit(1)
        
    if callset_csv.empty:
        # Generate an empty report if there are no calls, matching the standard format
        with open(args.output, 'w') as fh:
            fh.write("========================================================================\n")
            fh.write(f" CNV Caller Performance Evaluation Report - {args.workflow_type}\n")
            fh.write("========================================================================\n")
            fh.write(f"Source Callset CSV  : {args.input_csv}\n")
            fh.write(f"Truthset BED File   : {args.truth_bed}\n")
            fh.write(f"Probes Target BED   : {args.probes_bed}\n")
            fh.write(f"Sample Map File     : {args.map_file}\n")
            fh.write(f"Cohort Size         : {len(cohort_egan_samples)} Samples (mapped successfully)\n")
            fh.write("Called Rare CNVs    : 0 (inside this class)\n")
            fh.write("Truth set rare CNVs : 0 (within cohort)\n")
            fh.write("------------------------------------------------------------------------\n\n")
            fh.write("RIGOROUS SAMPLE-AWARE CNV CALLER EVALUATION METRICS\n")
            fh.write("------------------------------------------------------------------------\n")
            fh.write("  True Positives (TP)  : 0\n")
            fh.write("  False Positives (FP) : 0\n")
            fh.write("  False Negatives (FN) : 0\n")
            fh.write("  True Negatives (TN)  : 0\n\n")
            fh.write("  Sensitivity (Recall) : 0.000000 (0.00%)\n")
            fh.write("  Precision (PPV)      : 0.000000 (0.00%)\n")
            fh.write("  Specificity          : 0.000000 (0.00%)\n")
            fh.write("  F-beta (beta=2)      : 0.000000\n")
            fh.write("  Matthews Corr Coeff  : 0.000000\n")
            fh.write("========================================================================\n")
        logging.warning("Input CSV is empty. Wrote standard formatted empty report.")
        sys.exit(0)

    # Clean and map called CNVs to EGAN sample IDs (strip .md or other suffixes)
    callset_csv['SAMPLE_clean'] = callset_csv['SAMPLE'].astype(str).str.replace(r'\.(md|bam|sorted|normalised).*$', '', regex=True)
    callset_csv['sample'] = callset_csv['SAMPLE_clean'].map(egaf_to_egan)
    # Keep only rows that successfully mapped
    callset_mapped = callset_csv.dropna(subset=['sample']).copy()
    
    # Convert callset to 5-column BED representation
    callset_cnvs = callset_mapped[['CHROM', 'POS', 'END', 'SVTYPE', 'sample']].rename(
        columns={'CHROM': 'chr', 'POS': 'start', 'END': 'end', 'SVTYPE': 'cnv_type'}
    )
    callset_cnvs['start'] = callset_cnvs['start'].astype(int)
    callset_cnvs['end'] = callset_cnvs['end'].astype(int)
    callset_cnvs['chr'] = callset_cnvs['chr'].astype(str)
    callset_cnvs['sample'] = callset_cnvs['sample'].astype(str).str.strip()

    # 3. Load truth BED and filter for only our 90 cohort EGAN samples
    truth_all = load_truth_bed(args.truth_bed)
    truth_cnvs = truth_all[truth_all['sample'].isin(cohort_egan_samples)].copy()

    # 4. Load target probes BED
    probes = load_probes_bed(args.probes_bed)

    # ==========================================
    # EVALUATION MODE A: Rigorous Sample-Aware
    # ==========================================
    sample_tp = 0
    sample_fn = 0
    sample_fp = 0
    sample_tn = 0

    # Pre-group by sample for performance
    truth_by_sample = {s: grp for s, grp in truth_cnvs.groupby('sample')}
    callset_by_sample = {s: grp for s, grp in callset_cnvs.groupby('sample')}

    for s_egan in cohort_egan_samples:
        truth_s = truth_by_sample.get(s_egan, pd.DataFrame(columns=['chr', 'start', 'end']))
        callset_s = callset_by_sample.get(s_egan, pd.DataFrame(columns=['chr', 'start', 'end']))
        
        truth_overlaps = find_overlapping_probe_indices(probes, truth_s)
        callset_overlaps = find_overlapping_probe_indices(probes, callset_s)
        
        tp = len(truth_overlaps & callset_overlaps)
        fn = len(truth_overlaps - callset_overlaps)
        fp = len(callset_overlaps - truth_overlaps)
        tn = len(probes) - (tp + fn + fp)
        
        sample_tp += tp
        sample_fn += fn
        sample_fp += fp
        sample_tn += tn

    sens_s, prec_s, spec_s, f2_s, mcc_s = compute_metrics(sample_tp, sample_fn, sample_fp, sample_tn, beta=2)

    # 5. Write beautiful performance report
    with open(args.output, 'w') as fh:
        fh.write("========================================================================\n")
        fh.write(f" CNV Caller Performance Evaluation Report - {args.workflow_type}\n")
        fh.write("========================================================================\n")
        fh.write(f"Source Callset CSV  : {args.input_csv}\n")
        fh.write(f"Truthset BED File   : {args.truth_bed}\n")
        fh.write(f"Probes Target BED   : {args.probes_bed}\n")
        fh.write(f"Sample Map File     : {args.map_file}\n")
        fh.write(f"Cohort Size         : {len(cohort_egan_samples)} Samples (mapped successfully)\n")
        fh.write(f"Called Rare CNVs    : {len(callset_cnvs)} (inside this class)\n")
        fh.write(f"Truth set rare CNVs : {len(truth_cnvs)} (within cohort)\n")
        fh.write("------------------------------------------------------------------------\n\n")

        fh.write("RIGOROUS SAMPLE-AWARE CNV CALLER EVALUATION METRICS\n")
        fh.write("------------------------------------------------------------------------\n")
        fh.write("This evaluation assesses each sample individually at the probe level and\n")
        fh.write("sums the confusion matrix counts across the cohort, preventing false\n")
        fh.write("positives and false negatives from collapsing across different samples.\n\n")
        fh.write(f"  True Positives (TP)  : {sample_tp}\n")
        fh.write(f"  False Positives (FP) : {sample_fp}\n")
        fh.write(f"  False Negatives (FN) : {sample_fn}\n")
        fh.write(f"  True Negatives (TN)  : {sample_tn}\n\n")
        fh.write(f"  Sensitivity (Recall) : {sens_s:.6f} ({sens_s*100:.2f}%)\n")
        fh.write(f"  Precision (PPV)      : {prec_s:.6f} ({prec_s*100:.2f}%)\n")
        fh.write(f"  Specificity          : {spec_s:.6f} ({spec_s*100:.2f}%)\n")
        fh.write(f"  F-beta (beta=2)      : {f2_s:.6f}\n")
        fh.write(f"  Matthews Corr Coeff  : {mcc_s:.6f}\n")
        fh.write("========================================================================\n")

    logging.info(f"Successfully wrote performance report to {args.output}")

if __name__ == '__main__':
    main()
