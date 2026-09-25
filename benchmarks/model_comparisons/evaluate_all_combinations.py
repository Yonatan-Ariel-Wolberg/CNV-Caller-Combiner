#!/usr/bin/env python3
import os
import sys
import glob
import itertools
import argparse
import pandas as pd
from bisect import bisect_left

# Reuse loading/search functions from evaluate_caller_performance.py
def load_bed_file(bed_file):
    return pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end'],
                       usecols=[0, 1, 2])

def load_truth_bed(truth_file):
    # truth bed columns: chr, start, end, cnv_type, sample
    return pd.read_csv(truth_file, sep='\t', header=None,
                       names=['chr', 'start', 'end', 'cnv_type', 'sample'])

def parse_vcf_to_intervals(vcf_path):
    """
    Parses a single VCF/VCF.gz file and extracts (chr, start, end) intervals.
    """
    import gzip
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    intervals = []
    sample_id = None
    try:
        with open_func(vcf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    cols = line.strip().split('\t')
                    if len(cols) > 9:
                        sample_id = cols[9]
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 8:
                    continue
                # Skip reference calls (non-variant segments)
                alt = cols[4]
                if alt == '.' or alt == '':
                    continue
                chrom = cols[0]
                pos = int(cols[1])
                start = pos - 1
                info = cols[7]
                
                # Extract END coordinate
                end = start + 1
                for item in info.split(';'):
                    if item.startswith('END='):
                        try:
                            end = int(item.split('=')[1])
                        except ValueError:
                            pass
                intervals.append((chrom, start, end))
    except Exception as e:
        print(f"Warning: Failed to parse VCF {vcf_path}: {e}")
    
    if sample_id is None:
        # fallback to filename
        sample_id = os.path.basename(vcf_path).split('.')[0]
    return sample_id, intervals

def find_overlapped_probe_indices(probes_by_chrom, intervals):
    """
    Finds probe indices overlapping with a list of intervals.
    """
    overlapped = set()
    for chrom_raw, c_start, c_end in intervals:
        chrom = str(chrom_raw).replace('chr', '').replace('CHR', '').strip()
        if chrom not in probes_by_chrom:
            continue
        p_list, p_starts = probes_by_chrom[chrom]
        n_probes = len(p_list)
        
        idx = bisect_left(p_starts, c_start)
        
        # Check backwards
        for i in range(idx - 1, -1, -1):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                continue
            if p_end > c_start:
                overlapped.add(p_idx)
            if p_start < c_start - 50000:
                break
                
        # Check forwards
        for i in range(idx, n_probes):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                break
            if p_end > c_start:
                overlapped.add(p_idx)
                
    return overlapped

def main():
    parser = argparse.ArgumentParser(description="Evaluate all sub-combinations of CNV callers at probe level.")
    parser.add_argument("--probes_bed", required=True, help="Target probes BED file")
    parser.add_argument("--truth_bed", required=True, help="Truth set BED file")
    parser.add_argument("--sample_map", help="Sample mapping file")
    parser.add_argument("--results_dir", default="results", help="Nextflow results output directory")
    parser.add_argument("--output_report", default="results/comparative_combinations_report.md", help="Output markdown file path")
    args = parser.parse_args()
    
    # Load sample map if provided
    sample_map = {}
    if args.sample_map and os.path.exists(args.sample_map):
        print(f"Loading sample map from {args.sample_map}...")
        with open(args.sample_map, 'r') as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sample_map[parts[0]] = parts[1]
                    
    print("Loading probes and truth sets...")
    probes_df = load_bed_file(args.probes_bed)
    # Pre-structure probes for fast binary search
    probes_df['chr_str'] = probes_df['chr'].astype(str).str.replace(r'^[cC][hH][rR]', '', regex=True).str.strip()
    probes_by_chrom = {}
    for chrom, group in probes_df.groupby('chr_str'):
        sorted_probes = sorted([(row.Index, row.start, row.end) for row in group.itertuples()], key=lambda x: x[1])
        probes_by_chrom[chrom] = (sorted_probes, [p[1] for p in sorted_probes])
        
    truth_df = load_truth_bed(args.truth_bed)
    
    callers = {
        'CANOES': f"{args.results_dir}/out_CANOES/vcfs/*.vcf*",
        'CLAMMS': f"{args.results_dir}/out_CLAMMS/vcfs/*.vcf*",
        'XHMM': f"{args.results_dir}/out_XHMM/vcfs/*.vcf*",
        'CNVKIT': f"{args.results_dir}/out_CNVKIT/vcfs/*.vcf*",
        'GCNV': f"{args.results_dir}/out_GCNV/vcfs/*.vcf*",
        'DRAGEN': f"{args.results_dir}/out_DRAGEN/vcfs/*.vcf*"
    }
    
    # We will gather parsed data: sample -> caller -> set of called probe indices
    sample_caller_probes = {}
    
    print("Parsing VCFs for each caller...")
    all_samples = set()
    for caller_name, pattern in callers.items():
        vcf_files = glob.glob(pattern)
        if not vcf_files:
            # Try alternate paths or uncompressed/compressed files
            alt_pattern = pattern.replace('vcfs/*.vcf*', '*.vcf*')
            vcf_files = glob.glob(alt_pattern)
        
        print(f"  {caller_name}: found {len(vcf_files)} files")
        for vcf_path in vcf_files:
            # Skip indices (.tbi)
            if vcf_path.endswith('.tbi'):
                continue
            sample_id, intervals = parse_vcf_to_intervals(vcf_path)
            # Clean sample ID
            sample_id = sample_id.split('_')[0].split('.')[0]
            mapped_sample_id = sample_map.get(sample_id, sample_id)
            all_samples.add(mapped_sample_id)
            
            overlapped_probes = find_overlapped_probe_indices(probes_by_chrom, intervals)
            
            if mapped_sample_id not in sample_caller_probes:
                sample_caller_probes[mapped_sample_id] = {}
            sample_caller_probes[mapped_sample_id][caller_name] = overlapped_probes
            
    print(f"Total samples processed: {len(all_samples)}")
    
    # Map truth CNVs to target probes per sample
    print("Mapping truth set to target probes...")
    sample_truth_probes = {}
    for sample_id, group in truth_df.groupby('sample'):
        intervals = [(row.chr, row.start, row.end) for row in group.itertuples()]
        overlapped_probes = find_overlapped_probe_indices(probes_by_chrom, intervals)
        sample_truth_probes[sample_id] = overlapped_probes
        
    # Get all actual active callers
    active_callers = sorted(list(set(c for s in sample_caller_probes.values() for c in s.keys())))
    print(f"Active callers: {active_callers}")
    
    # Total targets list
    all_probe_indices = set(probes_df.index)
    
    # Evaluate combinations
    results = []
    
    # Combinations to test: subsets of callers from size 1 to N
    for r in range(1, len(active_callers) + 1):
        for combo in itertools.combinations(active_callers, r):
            combo_name = "+".join(combo)
            
            # 1. Union Merge (Called by >= 1)
            tp_union, fp_union, fn_union = 0, 0, 0
            # 2. Consensus Merge (Called by >= 2, only for combo size >= 2)
            tp_cons, fp_cons, fn_cons = 0, 0, 0
            
            for sample_id in all_samples:
                truth_set = sample_truth_probes.get(sample_id, set())
                
                # Get called probe sets for the selected callers for this sample
                caller_sets = [sample_caller_probes[sample_id].get(c, set()) for c in combo]
                
                # Union set
                union_set = set().union(*caller_sets)
                
                # Consensus set (support >= 2)
                consensus_set = set()
                if len(combo) >= 2:
                    # Count how many callers call each probe
                    counts = {}
                    for s in caller_sets:
                        for idx in s:
                            counts[idx] = counts.get(idx, 0) + 1
                    consensus_set = {idx for idx, count in counts.items() if count >= 2}
                
                # Evaluate Union
                tp_u = len(union_set.intersection(truth_set))
                fp_u = len(union_set.difference(truth_set))
                fn_u = len(truth_set.difference(union_set))
                tp_union += tp_u
                fp_union += fp_u
                fn_union += fn_u
                
                # Evaluate Consensus
                if len(combo) >= 2:
                    tp_c = len(consensus_set.intersection(truth_set))
                    fp_c = len(consensus_set.difference(truth_set))
                    fn_c = len(truth_set.difference(consensus_set))
                    tp_cons += tp_c
                    fp_cons += fp_c
                    fn_cons += fn_c
            
            # Helper to calculate stats
            def calc_metrics(tp, fp, fn):
                prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                f1 = 2 * prec * sens / (prec + sens) if (prec + sens) > 0 else 0.0
                return prec, sens, f1
                
            p_u, s_u, f_u = calc_metrics(tp_union, fp_union, fn_union)
            results.append({
                'Combination': combo_name,
                'Mode': 'Union (>=1)',
                'Size': len(combo),
                'TP': tp_union,
                'FP': fp_union,
                'FN': fn_union,
                'Precision': p_u,
                'Sensitivity': s_u,
                'F1-Score': f_u
            })
            
            if len(combo) >= 2:
                p_c, s_c, f_c = calc_metrics(tp_cons, fp_cons, fn_cons)
                results.append({
                    'Combination': combo_name,
                    'Mode': 'Consensus (>=2)',
                    'Size': len(combo),
                    'TP': tp_cons,
                    'FP': fp_cons,
                    'FN': fn_cons,
                    'Precision': p_c,
                    'Sensitivity': s_c,
                    'F1-Score': f_c
                })
                
    df_res = pd.DataFrame(results)
    df_res = df_res.sort_values(by='F1-Score', ascending=False)
    
    # Format Precision, Sensitivity, F1-Score as percentages
    df_res['Precision'] = (df_res['Precision'] * 100).round(2).astype(str) + '%'
    df_res['Sensitivity'] = (df_res['Sensitivity'] * 100).round(2).astype(str) + '%'
    df_res['F1-Score'] = (df_res['F1-Score'] * 100).round(2).astype(str) + '%'
    
    # Save markdown report
    os.makedirs(os.path.dirname(args.output_report), exist_ok=True)
    with open(args.output_report, 'w') as f:
        f.write("# Comprehensive Caller Combinations Performance Report\n\n")
        f.write("This report evaluates all possible combinations of individual CNV callers (union and consensus) at the probe level against the Decipher truth set.\n\n")
        f.write("## Performance Summary Table\n\n")
        
        headers = list(df_res.columns)
        f.write("| " + " | ".join(map(str, headers)) + " |\n")
        f.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
        for row in df_res.itertuples(index=False):
            f.write("| " + " | ".join(map(str, row)) + " |\n")
        f.write("\n")
        
    print(f"Report written to {args.output_report}")

if __name__ == '__main__':
    main()
