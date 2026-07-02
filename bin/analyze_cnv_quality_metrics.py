#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: analyze_cnv_quality_metrics.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import os
import argparse
import glob
import pandas as pd
import numpy as np

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
    records = []
    if not os.path.exists(vcf_path):
        return records
        
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
            
            # Parse INFO for END
            info = parts[7]
            end = pos
            for field in info.split(';'):
                if field.startswith('END='):
                    try:
                        end = int(field.split('=')[1])
                    except:
                        pass
            
            # Parse FORMAT and sample data
            fmt_keys = parts[8].split(':')
            fmt_vals = parts[9].split(':')
            fmt_map = dict(zip(fmt_keys, fmt_vals))
            
            # Extract metrics
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
                    
            records.append({
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
    return records

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
    if probes_df is None:
        return None
    from bisect import bisect_left
    probe_dict = {}
    for chrom, group in probes_df.groupby('CHROM_std'):
        probe_dict[chrom] = sorted(list(zip(group['START'], group['END'])))
        
    def count_overlapping_probes(chrom_std, pos, end):
        if chrom_std not in probe_dict:
            return 0
        chrom_probes = probe_dict[chrom_std]
        starts = [p[0] for p in chrom_probes]
        idx_end = bisect_left(starts, end)
        
        overlap_count = 0
        for i in range(idx_end):
            p_start, p_end = chrom_probes[i]
            if p_start < end and p_end > pos:
                overlap_count += 1
        return overlap_count
    return count_overlapping_probes

def main():
    args = parse_args()
    
    # 1. Load input CSV
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
    
    # 2. Load SVTYPE mapping from AnnotSV output
    svtype_map = None
    if os.path.exists(args.annotsv_tsv):
        try:
            annotsv_full = pd.read_csv(args.annotsv_tsv, sep='\t', keep_default_na=False, low_memory=False)
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
    
    # 3. Load probes and setup overlap counter
    probes_df = load_probes(args.probes_bed)
    probe_counter = build_probe_counter(probes_df)
    
    # 4. Extract metrics from individual VCF files
    vcf_cache = {}
    extracted_data = []
    
    for sample, sample_group in df.groupby('SAMPLE'):
        # Check in vcf_dir for standard filenames
        vcf_path = os.path.join(args.vcf_dir, f"{sample}.md_filtered.cnv.vcf")
        if not os.path.exists(vcf_path):
            vcf_path = os.path.join(args.vcf_dir, f"{sample}_filtered.cnv.vcf")
            
        if not os.path.exists(vcf_path):
            # Try VCF.GZ if present
            vcf_path_gz = vcf_path + ".gz"
            if os.path.exists(vcf_path_gz):
                vcf_path = vcf_path_gz
                
        if vcf_path not in vcf_cache:
            vcf_cache[vcf_path] = parse_vcf(vcf_path)
            
        vcf_records = vcf_cache[vcf_path]
        
        for idx, row in sample_group.iterrows():
            chrom = row['CHROM_std']
            pos = row['POS']
            end = row['END']
            
            # Match coordinate boundaries
            best_match = None
            best_overlap = -1
            for rec in vcf_records:
                if rec['chrom'] != chrom:
                    continue
                if abs(rec['pos'] - pos) <= 1000 and abs(rec['end'] - end) <= 1000:
                    overlap_start = max(rec['pos'], pos)
                    overlap_end = min(rec['end'], end)
                    overlap = max(0, overlap_end - overlap_start)
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best_match = rec
                        
            rec_data = {
                'CHROM': row['CHROM'],
                'POS': pos,
                'END': end,
                'SAMPLE': sample,
                'ACMG_class': row.get('ACMG_class', 'NA'),
                'SVTYPE': row['SVTYPE'],
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
    res_df.to_csv(args.output_csv, index=False)
    
    # 5. Compute CNVs per sample per ACMG class
    sample_acmg = res_df.groupby(['SAMPLE', 'ACMG_class']).size().unstack(fill_value=0)
    sample_acmg_csv = args.output_csv.replace("_cnv_vcf_metrics.csv", "_cnv_counts_per_sample_per_acmg.csv")
    sample_acmg.to_csv(sample_acmg_csv)
    
    # 6. Format Summary Report
    with open(args.output_summary, 'w') as out:
        out.write("========================================================================\n")
        out.write(f" CNV Caller Quality Metrics Assessment Report - {args.run_name}\n")
        out.write("========================================================================\n")
        out.write(f"Source Callset CSV  : {args.input_csv}\n")
        out.write(f"Target Probes BED   : {args.probes_bed}\n")
        out.write(f"Total Rare CNVs     : {len(res_df)}\n")
        out.write("------------------------------------------------------------------------\n\n")
        
        # CNV Type Distribution
        counts = res_df['SVTYPE'].value_counts()
        out.write("CNV TYPE DISTRIBUTION:\n")
        for svtype, cnt in counts.items():
            pct = (cnt / len(res_df)) * 100
            out.write(f"  {svtype:8s} : {cnt:5d} ({pct:.2f}%)\n")
        out.write("\n")
        
        # Size statistics
        out.write("SIZE STATISTICS BY CNV TYPE:\n")
        for svtype in res_df['SVTYPE'].unique():
            sub = res_df[res_df['SVTYPE'] == svtype]
            out.write(f"  {svtype}:\n")
            out.write(f"    Mean   : {sub['Size'].mean():,.1f} bp\n")
            out.write(f"    Median : {sub['Size'].median():,.1f} bp\n")
            out.write(f"    Range  : {sub['Size'].min():,} - {sub['Size'].max():,} bp\n")
        out.write("\n")
        
        # Probe overlaps
        out.write("PROBE OVERLAP COUNTS BY CNV TYPE:\n")
        for svtype in res_df['SVTYPE'].unique():
            sub = res_df[res_df['SVTYPE'] == svtype]
            out.write(f"  {svtype}:\n")
            out.write(f"    Mean   : {sub['Probe_Count'].mean():.2f} probes\n")
            out.write(f"    Median : {sub['Probe_Count'].median():.1f} probes\n")
            out.write(f"    Range  : {sub['Probe_Count'].min():.0f} - {sub['Probe_Count'].max():.0f} probes\n")
        out.write("\n")
        
        # CNVs per sample per ACMG class summary stats
        out.write("CNV CALLS PER SAMPLE BY ACMG CLASS (SUMMARY STATISTICS):\n")
        out.write(sample_acmg.describe().to_string() + "\n\n")
        
        # VCF Metrics statistics (describe summaries)
        out.write("VCF QUALITY FIELDS SUMMARY STATISTICS:\n")
        metrics = ['QUAL', 'BC', 'SM', 'LR', 'PE']
        metric_desc = res_df[metrics].describe().to_string()
        out.write(metric_desc + "\n\n")
        
        # Split metrics summaries by SVTYPE
        for metric in metrics:
            out.write(f"--- {metric} DISTRIBUTION BY SVTYPE ---\n")
            grouped = res_df.groupby('SVTYPE')[metric].describe().to_string()
            out.write(grouped + "\n\n")
            
        out.write("========================================================================\n")
        
    print(f"Metrics summary written to {args.output_summary}")

if __name__ == '__main__':
    main()
