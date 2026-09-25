#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: analyze_cnv_quality_metrics.py (OPTIMIZED)
# Purpose: This pipeline helper script extracts, compiles, and assesses CNV call-level
#   quality metrics from filtered VCF files for each sample, maps SVTYPE annotations
#   from AnnotSV, and counts overlapping probes from a target design BED file.
# =============================================================================

import os
import argparse
import glob
import pandas as pd
import numpy as np
from bisect import bisect_left

def parse_args():
    parser = argparse.ArgumentParser(description="Extract and assess CNV call-level VCF quality metrics.")
    parser.add_argument('--input-csv', required=True, help="Path to the annotated rare CNVs CSV file.")
    parser.add_argument('--annotsv-tsv', required=True, help="Path to the AnnotSV TSV file (for recovering SVTYPE).")
    parser.add_argument('--vcf-dir', required=True, help="Directory containing the individual sample filtered VCF files.")
    parser.add_argument('--probes-bed', required=True, help="Path to the target probes BED file.")
    parser.add_argument('--run-name', required=True, help="Name of the pipeline run (e.g. Control or Mosdepth).")
    parser.add_argument('--output-summary', required=True, help="Path to write the formatted summary text report.")
    parser.add_argument('--output-csv', required=True, help="Path to write the detailed call-level metrics CSV.")
    return parser.parse_args()

def standardize_chrom(c):
    c_str = str(c).strip()
    if c_str.lower().startswith('chr'):
        return c_str[3:]
    return c_str

def parse_vcf(vcf_path):
    """
    Parses a single sample CNV VCF file to extract records grouped by chromosome and sorted by position.
    
    Args:
        vcf_path (str): Path to the VCF file.
        
    Returns:
        dict: A dictionary mapping chrom_std -> list of dicts (sorted by pos).
    """
    records_by_chrom = {}
    if not os.path.exists(vcf_path):
        return records_by_chrom
        
    parsed_records = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
                
            chrom = standardize_chrom(parts[0])
            pos = int(parts[1])
            qual = float(parts[5]) if parts[5] != '.' else None
            
            info = parts[7]
            end = pos
            for field in info.split(';'):
                if field.startswith('END='):
                    try:
                        end = int(field.split('=')[1])
                    except:
                        pass
            
            fmt_keys = parts[8].split(':')
            fmt_vals = parts[9].split(':')
            fmt_map = dict(zip(fmt_keys, fmt_vals))
            
            sm = float(fmt_map['SM']) if 'SM' in fmt_map and fmt_map['SM'] != '.' else None
            cn = int(fmt_map['CN']) if 'CN' in fmt_map and fmt_map['CN'] != '.' else None
            bc = int(fmt_map['BC']) if 'BC' in fmt_map and fmt_map['BC'] != '.' else None
            lr = float(fmt_map['LR']) if 'LR' in fmt_map and fmt_map['LR'] != '.' else None
            
            pe_val = 0
            if 'PE' in fmt_map and fmt_map['PE'] != '.':
                try:
                    pe_val = sum(map(int, fmt_map['PE'].split(',')))
                except:
                    pass
                    
            parsed_records.append({
                'chrom': chrom,
                'pos': pos,
                'end': end,
                'qual': qual,
                'sm': sm,
                'cn': cn,
                'bc': bc,
                'lr': lr,
                'pe': pe_val
            })
            
    # Group by chromosome
    for rec in parsed_records:
        records_by_chrom.setdefault(rec['chrom'], []).append(rec)
        
    # Sort each chromosome list by start position
    for chrom in records_by_chrom:
        records_by_chrom[chrom].sort(key=lambda x: x['pos'])
        
    return records_by_chrom

def load_probes(probes_bed_path):
    if os.path.exists(probes_bed_path):
        try:
            probes_df = pd.read_csv(probes_bed_path, sep='\t', header=None, names=['CHROM', 'START', 'END'], usecols=[0, 1, 2])
            probes_df['CHROM_std'] = probes_df['CHROM'].apply(standardize_chrom)
            return probes_df
        except Exception as e:
            print(f"Error loading probes BED: {e}")
    return None

def build_probe_counter(probes_df):
    """
    Builds a high-performance bisect-based helper function to count probe overlaps.
    Precalculates interval lists and start indexes once during setup.
    """
    if probes_df is None:
        return None
    
    # Organize probes by standardized chromosome in sorted order for fast lookup
    probe_dict = {}
    probe_starts = {}
    for chrom, group in probes_df.groupby('CHROM_std'):
        sorted_probes = sorted(list(zip(group['START'], group['END'])))
        probe_dict[chrom] = sorted_probes
        probe_starts[chrom] = [p[0] for p in sorted_probes]
        
    def count_overlapping_probes(chrom_std, pos, end):
        if chrom_std not in probe_dict:
            return 0
        chrom_probes = probe_dict[chrom_std]
        starts = probe_starts[chrom_std]
        
        # Binary search for the first probe starting after or at the CNV end coordinate
        idx_end = bisect_left(starts, end)
        
        # Count probes before this index that overlap the CNV query window
        overlap_count = 0
        for i in range(idx_end - 1, -1, -1):
            p_start, p_end = chrom_probes[i]
            if p_start < end and p_end > pos:
                overlap_count += 1
            # Stop early if probe starts way before the CNV window starts (50kb target buffer)
            if p_start < pos - 50000:
                break
        return overlap_count
    return count_overlapping_probes

def main():
    args = parse_args()
    
    df = pd.read_csv(args.input_csv, keep_default_na=False)
    if df.empty:
        print("Input CSV is empty. Writing empty outputs.")
        with open(args.output_summary, 'w') as out:
            out.write("No CNVs detected for quality metrics analysis.\n")
        pd.DataFrame().to_csv(args.output_csv, index=False)
        return
        
    df['POS'] = df['POS'].astype(int)
    df['END'] = df['END'].astype(int)
    df['CHROM_std'] = df['CHROM'].apply(standardize_chrom)
    df['Size'] = df['END'] - df['POS']
    
    svtype_map = None
    if os.path.exists(args.annotsv-tsv if hasattr(args, 'annotsv-tsv') else args.annotsv_tsv):
        annotsv_path = args.annotsv_tsv
        try:
            annotsv_full = pd.read_csv(annotsv_path, sep='\t', keep_default_na=False, low_memory=False)
            annotsv_full = annotsv_full[annotsv_full['Annotation_mode'] == 'full']
            annotsv_map = annotsv_full[['SV_chrom', 'SV_start', 'SV_end', 'SV_type']].copy()
            annotsv_map['CHROM_std'] = annotsv_map['SV_chrom'].apply(standardize_chrom)
            annotsv_map['POS'] = annotsv_map['SV_start'].astype(int)
            annotsv_map['END'] = annotsv_map['SV_end'].astype(int)
            annotsv_map = annotsv_map.rename(columns={'SV_type': 'AnnotSV_SV_type'})
            svtype_map = annotsv_map[['CHROM_std', 'POS', 'END', 'AnnotSV_SV_type']].drop_duplicates()
        except Exception as e:
            print(f"Error loading AnnotSV map: {e}")
            
    if svtype_map is not None:
        df = pd.merge(df, svtype_map, on=['CHROM_std', 'POS', 'END'], how='left')
        df['SVTYPE'] = df['AnnotSV_SV_type'].fillna(df['SVTYPE'])
        df = df.drop(columns=['AnnotSV_SV_type'])
        
    df['SVTYPE'] = df['SVTYPE'].replace({'NA': 'UNKNOWN', '': 'UNKNOWN'})
    
    probes_df = load_probes(args.probes_bed)
    probe_counter = build_probe_counter(probes_df)
    
    vcf_cache = {}
    extracted_data = []
    
    # Process each sample using fast itertuples loops
    for sample, sample_group in df.groupby('SAMPLE'):
        vcf_path = os.path.join(args.vcf_dir, f"{sample}.md_filtered.cnv.vcf")
        if not os.path.exists(vcf_path):
            vcf_path = os.path.join(args.vcf_dir, f"{sample}_filtered.cnv.vcf")
            
        if not os.path.exists(vcf_path):
            vcf_path_gz = vcf_path + ".gz"
            if os.path.exists(vcf_path_gz):
                vcf_path = vcf_path_gz
                
        if vcf_path not in vcf_cache:
            vcf_cache[vcf_path] = parse_vcf(vcf_path)
            
        vcf_records = vcf_cache[vcf_path]
        
        # Itertuples (10x faster than iterrows)
        for row in sample_group.itertuples(index=False):
            chrom = row.CHROM_std
            pos = row.POS
            end = row.END
            
            best_match = None
            best_overlap = -1
            
            # Retrieve chromosome-specific sorted list
            chrom_recs = vcf_records.get(chrom, [])
            if chrom_recs:
                # Binary search for the first VCF record starting within the buffer window
                starts = [r['pos'] for r in chrom_recs]
                idx = bisect_left(starts, pos - 1000)
                
                for i in range(idx, len(chrom_recs)):
                    rec = chrom_recs[i]
                    if rec['pos'] > pos + 1000:
                        break
                    if abs(rec['end'] - end) <= 1000:
                        overlap_start = max(rec['pos'], pos)
                        overlap_end = min(rec['end'], end)
                        overlap = max(0, overlap_end - overlap_start)
                        if overlap > best_overlap:
                            best_overlap = overlap
                            best_match = rec
                        
            rec_data = {
                'CHROM': row.CHROM,
                'POS': pos,
                'END': end,
                'SAMPLE': sample,
                'ACMG_class': getattr(row, 'ACMG_class', 'NA'),
                'SVTYPE': row.SVTYPE,
                'Size': end - pos,
                'QUAL': None,
                'SM': None,
                'CN': None,
                'BC': None,
                'LR': None,
                'PE': None,
                'Probe_Count': 0
            }
            
            if best_match:
                rec_data['QUAL'] = best_match['qual']
                rec_data['SM'] = best_match['sm']
                rec_data['CN'] = best_match['cn']
                rec_data['BC'] = best_match['bc']
                rec_data['LR'] = best_match['lr']
                rec_data['PE'] = best_match['pe']
                
            if probe_counter:
                rec_data['Probe_Count'] = probe_counter(chrom, pos, end)
                
            extracted_data.append(rec_data)
            
    res_df = pd.DataFrame(extracted_data)
    
    if res_df.empty:
        # Create empty CSVs and empty summary report to guarantee pipeline resilience
        empty_cols = [
            'SAMPLE', 'CNV_ID', 'CHROM', 'POS', 'END', 'SVTYPE', 'Size',
            'ACMG_class', 'QUAL', 'SM', 'CN', 'BC', 'LR', 'PE', 'Probe_Count'
        ]
        res_df = pd.DataFrame(columns=empty_cols)
        res_df.to_csv(args.output_csv, index=False)
        sample_acmg_csv = args.output_csv.replace("_cnv_vcf_metrics.csv", "_cnv_counts_per_sample_per_acmg.csv")
        pd.DataFrame().to_csv(sample_acmg_csv)

        with open(args.output_summary, 'w') as out:
            out.write("========================================================================\n")
            out.write(f" CNV Caller Quality Metrics Assessment Report - {args.run_name}\n")
            out.write("========================================================================\n")
            out.write(f"Source Callset CSV  : {args.input_csv}\n")
            out.write(f"Target Probes BED   : {args.probes_bed}\n")
            out.write(f"Total Rare CNVs     : 0\n")
            out.write("------------------------------------------------------------------------\n\n")
            out.write("No CNVs found in the input callset.\n")
            out.write("========================================================================\n")
        print(f"Metrics summary written to {args.output_summary}")
        return

    res_df.to_csv(args.output_csv, index=False)
    
    sample_acmg = res_df.groupby(['SAMPLE', 'ACMG_class']).size().unstack(fill_value=0)
    sample_acmg_csv = args.output_csv.replace("_cnv_vcf_metrics.csv", "_cnv_counts_per_sample_per_acmg.csv")
    sample_acmg.to_csv(sample_acmg_csv)
    
    with open(args.output_summary, 'w') as out:
        out.write("========================================================================\n")
        out.write(f" CNV Caller Quality Metrics Assessment Report - {args.run_name}\n")
        out.write("========================================================================\n")
        out.write(f"Source Callset CSV  : {args.input_csv}\n")
        out.write(f"Target Probes BED   : {args.probes_bed}\n")
        out.write(f"Total Rare CNVs     : {len(res_df)}\n")
        out.write("------------------------------------------------------------------------\n\n")
        
        counts = res_df['SVTYPE'].value_counts()
        out.write("CNV TYPE DISTRIBUTION:\n")
        for svtype, cnt in counts.items():
            pct = (cnt / len(res_df)) * 100
            out.write(f"  {svtype:8s} : {cnt:5d} ({pct:.2f}%)\n")
        out.write("\n")
        
        out.write("SIZE STATISTICS BY CNV TYPE:\n")
        for svtype in res_df['SVTYPE'].unique():
            sub = res_df[res_df['SVTYPE'] == svtype]
            out.write(f"  {svtype}:\n")
            out.write(f"    Mean   : {sub['Size'].mean():,.1f} bp\n")
            out.write(f"    Median : {sub['Size'].median():,.1f} bp\n")
            out.write(f"    Range  : {sub['Size'].min():,} - {sub['Size'].max():,} bp\n")
        out.write("\n")
        
        out.write("PROBE OVERLAP COUNTS BY CNV TYPE:\n")
        for svtype in res_df['SVTYPE'].unique():
            sub = res_df[res_df['SVTYPE'] == svtype]
            out.write(f"  {svtype}:\n")
            out.write(f"    Mean   : {sub['Probe_Count'].mean():.2f} probes\n")
            out.write(f"    Median : {sub['Probe_Count'].median():.1f} probes\n")
            out.write(f"    Range  : {sub['Probe_Count'].min():.0f} - {sub['Probe_Count'].max():.0f} probes\n")
        out.write("\n")
        
        out.write("CNV CALLS PER SAMPLE BY ACMG CLASS (SUMMARY STATISTICS):\n")
        out.write(sample_acmg.describe().to_string() + "\n\n")
        
        out.write("VCF QUALITY FIELDS SUMMARY STATISTICS:\n")
        metrics = ['QUAL', 'BC', 'SM', 'LR', 'PE']
        metric_desc = res_df[metrics].describe().to_string()
        out.write(metric_desc + "\n\n")
        
        for metric in metrics:
            out.write(f"--- {metric} DISTRIBUTION BY SVTYPE ---\n")
            grouped = res_df.groupby('SVTYPE')[metric].describe().to_string()
            out.write(grouped + "\n\n")
            
        out.write("========================================================================\n")
        
    print(f"Metrics summary written to {args.output_summary}")

if __name__ == '__main__':
    main()
