#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: filter-rare-cnvs-ddduk-90.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import argparse
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate internal AF for cohort CNVs and filter rare variants.")
    parser.add_argument('--input-tsv', required=True, help="Path to the flattened TSV from bcftools.")
    parser.add_argument('--output-csv', required=True, help="Path to write the final filtered CSV report.")
    parser.add_argument('--total-samples', type=int, required=True, help="Total number of samples in the cohort (e.g., 90).")
    parser.add_argument('--max-af', type=float, default=0.05, help="Maximum allowed internal allele frequency (default: 0.05 / 5%%).")
    parser.add_argument('--annotsv-tsv', help="Path to the AnnotSV TSV annotation report (optional).")
    parser.add_argument('--output-regions', help="Path to write the unique rare CNV regions (tab-separated CHROM, POS, END) (optional).")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Calculate the max number of samples that can have this CNV to pass the filter
    # Ensure at least singletons (count = 1) are kept in small cohorts
    max_allowed_count = max(1.0, args.total_samples * args.max_af)
    logging.info(f"Cohort size: {args.total_samples}. Max frequency: {args.max_af*100}%.")
    logging.info(f"A CNV must appear in <= {max_allowed_count} samples to be retained.")

    # 1. Load the TSV
    try:
        df = pd.read_csv(args.input_tsv, sep='\t', names=['CHROM', 'POS', 'END', 'SVTYPE', 'SAMPLE', 'GT'], keep_default_na=False)
    except Exception as e:
        logging.error(f"Failed to read input TSV: {e}")
        sys.exit(1)

    if df.empty:
        logging.warning("Input TSV is empty. No CNVs passed initial filtering.")
        # Create empty output files to prevent Nextflow missing output crashes
        pd.DataFrame(columns=['CHROM', 'POS', 'END', 'SVTYPE', 'COUNT', 'INTERNAL_AF', 'SAMPLE', 'GT']).to_csv(args.output_csv, index=False)
        if args.output_regions:
            open(args.output_regions, 'w').close()
        sys.exit(0)

    # 2. Group by the variant coordinates to count occurrences
    cnv_counts = df.groupby(['CHROM', 'POS', 'END', 'SVTYPE']).size().reset_index(name='COUNT')

    # 3. Calculate the Internal Allele Frequency
    cnv_counts['INTERNAL_AF'] = cnv_counts['COUNT'] / args.total_samples

    # 4. Filter for Rare CNVs
    rare_cnvs = cnv_counts[cnv_counts['COUNT'] <= max_allowed_count]
    
    logging.info(f"Total unique CNVs before filter: {len(cnv_counts)}")
    logging.info(f"Rare unique CNVs after filter: {len(rare_cnvs)}")

    # If output-regions is requested, write unique rare CNV coordinates
    if args.output_regions:
        try:
            regions_df = rare_cnvs[['CHROM', 'POS', 'END']].copy()
            regions_df['CHROM'] = regions_df['CHROM'].astype(str)
            regions_df = regions_df.sort_values(by=['CHROM', 'POS'])
            regions_df.to_csv(args.output_regions, sep='\t', index=False, header=False)
            logging.info(f"Successfully wrote rare CNV regions to {args.output_regions}")
        except Exception as e:
            logging.error(f"Failed to write rare regions file: {e}")

    # 5. Merge the sample names back onto the filtered CNVs
    final_report = pd.merge(rare_cnvs, df, on=['CHROM', 'POS', 'END', 'SVTYPE'])

    # Optional: Merge AnnotSV annotation
    if args.annotsv_tsv:
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
            annotsv_filtered['SVTYPE'] = annotsv_filtered['SVTYPE'].astype(str)

            # Create standardized CHROM key for merging (stripping 'chr' prefix)
            def standardize_chrom(c):
                c_str = str(c).strip()
                if c_str.lower().startswith('chr'):
                    return c_str[3:]
                return c_str

            final_report['CHROM_std'] = final_report['CHROM'].apply(standardize_chrom)
            annotsv_filtered['CHROM_std'] = annotsv_filtered['SV_chrom_std'].apply(standardize_chrom)

            # Merge
            final_report = pd.merge(
                final_report,
                annotsv_filtered,
                on=['CHROM_std', 'POS', 'END', 'SVTYPE'],
                how='left'
            )
            # Drop temporary merge keys
            final_report = final_report.drop(columns=['CHROM_std', 'SV_chrom_std'], errors='ignore')
            logging.info("Successfully merged AnnotSV annotations into final report.")
        except Exception as e:
            logging.error(f"Failed to load or merge AnnotSV annotations: {e}")

    # 6. Sort for readability (by Chromosome, then Position)
    # Ensure CHROM is treated as string for sorting (e.g., chr1, chr2, chrX)
    final_report['CHROM'] = final_report['CHROM'].astype(str) 
    final_report = final_report.sort_values(by=['CHROM', 'POS'])

    # 7. Save to CSV
    final_report.to_csv(args.output_csv, index=False)
    logging.info(f"Successfully wrote rare CNV report to {args.output_csv}")

if __name__ == '__main__':
    main()
