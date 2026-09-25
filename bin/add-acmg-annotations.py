#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: add-acmg-annotations.py
# Purpose: Merge AnnotSV ACMG annotations back into the rare CNVs report.
#   This pipeline helper script matches chromosomes and coordinate ranges
#   between a sample/cohort filtered rare CNVs CSV and an AnnotSV TSV output.
#   It then splits the annotated results into separate CSV files based on
#   ACMG classes (e.g., class 5, 4, 3, 2, 1, NA) and high-priority (classes 4 and 5).
# Inputs:
#   - --input-csv: Path to the filtered rare CNVs CSV file.
#   - --annotsv-tsv: Path to the AnnotSV TSV annotation report.
#   - --output-csv: Path to write the final ACMG annotated CSV report.
# Outputs:
#   - An annotated CSV report containing ACMG classification and AnnotSV metrics.
#   - Split CSV reports for individual ACMG classes: _ACMG5.csv, _ACMG4.csv, etc.
#   - A combined high-priority CSV report: _ACMG4_5.csv.
# Operational notes:
#   - Handled with fallback handling for empty inputs to maintain Nextflow compatibility.
#   - Ensures coordinate comparisons are performed after parsing fields as integers.
#   - Chromosome merging is standardized by stripping the 'chr' prefix.
# =============================================================================
import argparse
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    """
    Parses command-line arguments for merging AnnotSV ACMG annotations.
    
    Returns:
        argparse.Namespace: Parsed CLI arguments containing input/output paths.
    """
    parser = argparse.ArgumentParser(description="Merge AnnotSV ACMG annotations back into the rare CNVs report.")
    parser.add_argument('--input-csv', required=True, help="Path to the filtered rare CNVs CSV.")
    parser.add_argument('--annotsv-tsv', required=True, help="Path to the AnnotSV TSV annotation report.")
    parser.add_argument('--output-csv', required=True, help="Path to write the final ACMG annotated CSV report.")
    return parser.parse_args()

def main():
    """
    Main function to load input CSV, parse AnnotSV TSV, standardize coordinates/chromosomes,
    merge ACMG annotations, and split the final reports into classes.
    """
    args = parse_args()
    
    # 1. Load the rare CNVs CSV
    try:
        final_report = pd.read_csv(args.input_csv, keep_default_na=False)
        # Drop pre-existing ACMG columns to make the run idempotent
        cols_to_drop = ['ACMG_class', 'AnnotSV_ranking_score', 'AnnotSV_ranking_criteria', 'OVERLAPPED_GENES']
        final_report = final_report.drop(columns=cols_to_drop, errors='ignore')
    except Exception as e:
        logging.error(f"Failed to read input CSV: {e}")
        sys.exit(1)

    if final_report.empty:
        logging.warning("Input CSV is empty. No variants to annotate.")
        final_report.to_csv(args.output_csv, index=False)
        # Also write at least one split file to satisfy Nextflow glob outputs
        if args.output_csv.endswith('.csv'):
            out_path = args.output_csv[:-4] + "_ACMGNA.csv"
        else:
            out_path = args.output_csv + "_ACMGNA"
        final_report.to_csv(out_path, index=False)
        sys.exit(0)

    # 2. Load the AnnotSV TSV
    logging.info(f"Loading AnnotSV annotations from {args.annotsv_tsv}...")
    try:
        annotsv_df = pd.read_csv(args.annotsv_tsv, sep='\t', keep_default_na=False, low_memory=False)
        # Filter for full annotation lines
        annotsv_full = annotsv_df[annotsv_df['Annotation_mode'] == 'full'].copy()
        logging.info(f"Loaded {len(annotsv_full)} full annotations from AnnotSV.")

        # Keep only required columns
        cols_to_keep = {
            'SV_chrom': 'SV_chrom_std',
            'SV_start': 'POS',
            'SV_end': 'END',
            'SV_type': 'SVTYPE',
            'ACMG_class': 'ACMG_class',
            'AnnotSV_ranking_score': 'AnnotSV_ranking_score',
            'AnnotSV_ranking_criteria': 'AnnotSV_ranking_criteria',
            'Gene_name': 'OVERLAPPED_GENES'
        }
        # Intersect with available columns in AnnotSV output
        cols_to_keep = {k: v for k, v in cols_to_keep.items() if k in annotsv_full.columns}
        
        annotsv_filtered = annotsv_full[list(cols_to_keep.keys())].rename(columns=cols_to_keep)
        
        # Standardize coordinates
        annotsv_filtered['POS'] = annotsv_filtered['POS'].astype(int)
        annotsv_filtered['END'] = annotsv_filtered['END'].astype(int)
        if 'SVTYPE' in annotsv_filtered.columns:
            annotsv_filtered['SVTYPE'] = annotsv_filtered['SVTYPE'].astype(str)

        # Create standardized CHROM key for merging (stripping 'chr' prefix)
        def standardize_chrom(c):
            c_str = str(c).strip()
            if c_str.lower().startswith('chr'):
                return c_str[3:]
            return c_str

        final_report['CHROM_std'] = final_report['CHROM'].apply(standardize_chrom)
        annotsv_filtered['CHROM_std'] = annotsv_filtered['SV_chrom_std'].apply(standardize_chrom)

        # Drop SVTYPE from AnnotSV dataframe to avoid column collision, and merge on coordinates only
        # This is robust against cases where SVTYPE is missing/NA in the rare CSV
        if 'SVTYPE' in annotsv_filtered.columns:
            annotsv_filtered = annotsv_filtered.drop(columns=['SVTYPE'])

        # Deduplicate AnnotSV annotations by coordinates to avoid row multiplication
        annotsv_filtered = annotsv_filtered.drop_duplicates(subset=['CHROM_std', 'POS', 'END'])

        # Merge
        final_report = pd.merge(
            final_report,
            annotsv_filtered,
            on=['CHROM_std', 'POS', 'END'],
            how='left'
        )
        # Drop temporary merge keys
        final_report = final_report.drop(columns=['CHROM_std', 'SV_chrom_std'], errors='ignore')
        logging.info("Successfully merged AnnotSV annotations into final report.")
    except Exception as e:
        logging.error(f"Failed to load or merge AnnotSV annotations: {e}")

    # 3. Sort for readability (by Chromosome, then Position)
    final_report['CHROM'] = final_report['CHROM'].astype(str) 
    final_report = final_report.sort_values(by=['CHROM', 'POS'])

    # 4. Save to CSV
    final_report.to_csv(args.output_csv, index=False)
    logging.info(f"Successfully wrote ACMG annotated rare CNV report to {args.output_csv}")

    # 5. Split by ACMG class into separate files depending on what scores are present
    def get_acmg_group(val):
        if pd.isna(val):
            return 'NA'
        s = str(val).strip()
        if s == '' or s.lower() in ['na', 'nan', '<na>']:
            return 'NA'
        try:
            f = float(s)
            if f.is_integer():
                return str(int(f))
            return str(f)
        except ValueError:
            return s

    if 'ACMG_class' not in final_report.columns:
        final_report['ACMG_class'] = 'NA'

    grouped_acmg = final_report['ACMG_class'].apply(get_acmg_group)

    for acmg_val in ['5', '4', '3', '2', '1', 'NA']:
        sub_df = final_report[grouped_acmg == acmg_val]
        if args.output_csv.endswith('.csv'):
            out_path = args.output_csv[:-4] + f"_ACMG{acmg_val}.csv"
        else:
            out_path = args.output_csv + f"_ACMG{acmg_val}"
        sub_df.to_csv(out_path, index=False)
        logging.info(f"Wrote {len(sub_df)} CNVs with ACMG={acmg_val} to {out_path}")

    # 6. Save a separate file containing CNVs with ACMG=5 and ACMG=4
    high_priority_df = final_report[grouped_acmg.isin(['4', '5'])]
    if args.output_csv.endswith('.csv'):
        out_path_high = args.output_csv[:-4] + "_ACMG4_5.csv"
    else:
        out_path_high = args.output_csv + "_ACMG4_5"
    high_priority_df.to_csv(out_path_high, index=False)
    logging.info(f"Wrote {len(high_priority_df)} CNVs with ACMG=4 and ACMG=5 to {out_path_high}")


if __name__ == '__main__':
    main()
