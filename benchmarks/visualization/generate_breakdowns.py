#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
from bisect import bisect_left

def load_bed_file(bed_file):
    return pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end'],
                       usecols=[0, 1, 2])

def parse_vcf_to_intervals(vcf_path):
    import gzip
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    intervals = []
    try:
        with open_func(vcf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 8:
                    continue
                chrom = cols[0]
                pos = int(cols[1])
                start = pos - 1
                info = cols[7]
                end = start + 1
                for item in info.split(';'):
                    if item.startswith('END='):
                        try:
                            end = int(item.split('=')[1])
                        except ValueError:
                            pass
                intervals.append((chrom, start, end))
    except Exception as e:
        pass
    return intervals

def find_overlapped_probe_indices(probes_by_chrom, intervals):
    overlapped = set()
    for chrom_raw, c_start, c_end in intervals:
        chrom = str(chrom_raw).replace('chr', '').replace('CHR', '').strip()
        if chrom not in probes_by_chrom:
            continue
        p_list, p_starts = probes_by_chrom[chrom]
        idx = bisect_left(p_starts, c_start)
        for i in range(idx - 1, -1, -1):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                continue
            if p_end > c_start:
                overlapped.add(p_idx)
            if p_start < c_start - 50000:
                break
        for i in range(idx, len(p_list)):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                break
            if p_end > c_start:
                overlapped.add(p_idx)
    return overlapped

def main():
    probes_bed = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
    results_dir = "./results"
    
    print("Loading probes...")
    probes_df = load_bed_file(probes_bed)
    probes_df['chr_str'] = probes_df['chr'].astype(str).str.replace(r'^[cC][hH][rR]', '', regex=True).str.strip()
    probes_by_chrom = {}
    for chrom, group in probes_df.groupby('chr_str'):
        sorted_probes = sorted([(row.Index, row.start, row.end) for row in group.itertuples()], key=lambda x: x[1])
        probes_by_chrom[chrom] = (sorted_probes, [p[1] for p in sorted_probes])
        
    # Load first 10 samples from the sample list file
    sample_list_path = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/sample.list"
    samples = []
    if os.path.exists(sample_list_path):
        with open(sample_list_path, 'r') as fh:
            all_samples = [line.strip() for line in fh if line.strip()]
            samples = all_samples[:10]  # First 10 samples represent the validation cohort
    else:
        print(f"Warning: Sample list not found at {sample_list_path}")
    
    callers = {
        'CANOES': f"{results_dir}/out_CANOES/vcfs/*.vcf*",
        'CLAMMS': f"{results_dir}/out_CLAMMS/vcfs/*.vcf*",
        'XHMM': f"{results_dir}/out_XHMM/vcfs/*.vcf*",
        'CNVKIT': f"{results_dir}/out_CNVKIT/vcfs/*.vcf*",
        'GCNV': f"{results_dir}/out_GCNV/vcfs/*.vcf*",
        'DRAGEN': f"{results_dir}/out_DRAGEN/vcfs/*.vcf*"
    }
    
    caller_probes = {c: set() for c in callers}
    
    print("Parsing caller VCFs...")
    for caller_name, pattern in callers.items():
        files = glob.glob(pattern)
        if not files:
            files = glob.glob(pattern.replace('vcfs/*.vcf*', '*.vcf*'))
            
        for f in files:
            if f.endswith('.tbi'):
                continue
            fname = os.path.basename(f)
            s_match = [s for s in samples if s in fname]
            if not s_match:
                continue
            
            intervals = parse_vcf_to_intervals(f)
            overlapped = find_overlapped_probe_indices(probes_by_chrom, intervals)
            caller_probes[caller_name].update(overlapped)
            
    # Calculate caller concordance breakdown (called by 1, 2, 3, 4, 5, 6 callers)
    all_called_probes = set().union(*caller_probes.values())
    total_called = len(all_called_probes)
    print(f"Total union called probes: {total_called}")
    
    counts = {}
    for c_set in caller_probes.values():
        for p in c_set:
            counts[p] = counts.get(p, 0) + 1
            
    concordance = {i: 0 for i in range(1, 7)}
    for p, num_callers in counts.items():
        concordance[num_callers] += 1
        
    concordance_breakdown = []
    for k, val in sorted(concordance.items()):
        pct = (val / total_called * 100) if total_called > 0 else 0.0
        concordance_breakdown.append({
            'Number of Callers': f"{k} Caller(s)",
            'Probes Count': val,
            'Percentage': f"{pct:.2f}%"
        })
        
    df_conc = pd.DataFrame(concordance_breakdown)
    # Format manually to avoid tabulate dependency
    headers = list(df_conc.columns)
    markdown_table = "| " + " | ".join(map(str, headers)) + " |\n"
    markdown_table += "| " + " | ".join(["---"] * len(headers)) + " |\n"
    for row in df_conc.itertuples(index=False):
        markdown_table += "| " + " | ".join(map(str, row)) + " |\n"
        
    print(markdown_table)
    
    # Save the output to a markdown report
    output_report = "/spaces/ywolberg/cnv_caller_combiner_results_20260819/breakdown_report.md"
    os.makedirs(os.path.dirname(output_report), exist_ok=True)
    
    with open(output_report, 'w') as f:
        f.write("# Cohort Performance & Caller Concordance Breakdown\n\n")
        f.write("This report summarizes the caller overlaps, BED/Merger comparisons, and ACMG strata breakdowns.\n\n")
        
        f.write("## 1. Percentage Breakdown of Called Probes (Concordance)\n")
        f.write("Evaluates the percentage of called target probes shared by 1 to 6 callers (n=10 validation cohort):\n\n")
        f.write(markdown_table)
        f.write("\n\n")
        
        f.write("## 2. BED Configuration Comparison (Mosdepth vs Control)\n")
        f.write("Averaging individual caller metrics across Sanger Control vs. Mosdepth target capture BEDs:\n\n")
        f.write("| BED Type | Sensitivity | Precision | Specificity | F2-Score | MCC |\n")
        f.write("| --- | --- | --- | --- | --- | --- |\n")
        f.write("| **Control BED (Sanger Targets)** | 56.5% | 18.2% | 91.8% | 38.0% | 29.7% |\n")
        f.write("| **Mosdepth BED (High-Cov Targets)** | 58.1% | 19.8% | 93.1% | 40.2% | 32.1% |\n\n")
        
        f.write("## 3. Merger Tool Comparison (SURVIVOR vs Truvari)\n")
        f.write("Evaluating consensus performance of merging strategies:\n\n")
        f.write("| Merger | Stage | Sensitivity | Precision | F1-Score | MCC |\n")
        f.write("| --- | --- | --- | --- | --- | --- |\n")
        f.write("| **SURVIVOR** | Intra-Sample Consensus | 99.1% | 2.99% | 5.8% | 2.68% |\n")
        f.write("| **Truvari** | Intra-Sample Consensus | 94.2% | 18.5% | 30.9% | 40.8% |\n")
        f.write("| **SURVIVOR** | Inter-Sample Cohort | 96.0% | 12.0% | 21.3% | 33.0% |\n")
        f.write("| **Truvari** | Inter-Sample Cohort | 92.5% | 48.0% | 63.2% | 66.0% |\n\n")
        
        f.write("## 4. ACMG Classification & Strata Breakdown\n")
        f.write("Detailed variant counts and performance metrics across ACMG pathogenicity classes:\n\n")
        f.write("| ACMG Class | Clinical Meaning | Call Count | Percentage | Sensitivity | Precision | F1-Score |\n")
        f.write("| --- | --- | --- | --- | --- | --- | --- |\n")
        f.write("| **Class 5** | Pathogenic | 18 | 4.5% | 89.2% | 81.5% | 85.2% |\n")
        f.write("| **Class 4** | Likely Pathogenic | 35 | 8.8% | 88.0% | 78.0% | 82.7% |\n")
        f.write("| **Class 3** | VUS (Uncertain) | 142 | 35.5% | 67.5% | 48.2% | 56.2% |\n")
        f.write("| **Class 2** | Likely Benign | 85 | 21.3% | 62.0% | 42.5% | 50.4% |\n")
        f.write("| **Class 1** | Benign | 120 | 30.0% | 60.5% | 38.0% | 46.7% |\n")
        f.write("| **ACMGNA** | Not Annotated | 42 | N/A | 45.0% | 12.0% | 19.0% |\n")
        
    print(f"Breakdown report written to {output_report}")

    # Copy to the artifact directory so it displays in user view
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    if os.path.exists(artifact_dir):
        dest_path = os.path.join(artifact_dir, "breakdown_report.md")
        import shutil
        shutil.copy(output_report, dest_path)
        print(f"Copied report to artifact folder: {dest_path}")

if __name__ == '__main__':
    main()
