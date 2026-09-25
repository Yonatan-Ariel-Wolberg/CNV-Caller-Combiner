#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: filter-rare-cnvs-ddduk-90.py
# Purpose: Calculate internal allele frequency (AF) for cohort CNVs and filter
#   out common variants to output rare CNVs (frequency <= max-af).
#   IMPORTANT: Allele frequency filtering focuses strictly and exclusively on the
#   internal dataset/cohort samples. It does NOT use external databases
#   (such as ClinVar, gnomAD, 1000 Genomes, or DGV).
# Inputs:
#   - --input-tsv: Flattened TSV file from bcftools (CHROM, POS, END, SVTYPE, SAMPLE, GT).
#   - --output-csv: Path to write the final filtered CSV report.
#   - --total-samples: Cohort size (integer). Used to calculate internal frequency.
#   - --max-af: Maximum internal allele frequency threshold (float, default: 0.05).
#   - --annotsv-tsv: Deprecated/ignored; external databases are not used for AF filtering.
#   - --output-regions: Path to write unique rare CNV coordinate regions (optional).
# Outputs:
#   - Filtered CSV containing rare CNVs with calculated internal COUNT and INTERNAL_AF metrics.
#   - Optional regions file containing tab-separated CHROM, POS, END columns of rare CNVs.
# Operational notes:
#   - At least singletons (count = 1) are kept in small cohorts even if they exceed max-af.
#   - Reciprocal overlap (>= 50%) is used to cluster overlapping calls across samples.
# =============================================================================
import argparse
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    """
    Parses command-line arguments for cohort CNV frequency filtering.
    
    Returns:
        argparse.Namespace: CLI parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate internal AF for cohort CNVs and filter rare variants (internal dataset focus only)."
    )
    parser.add_argument('--input-tsv', required=True, help="Path to the flattened TSV from bcftools.")
    parser.add_argument('--output-csv', required=True, help="Path to write the final filtered CSV report.")
    parser.add_argument('--total-samples', type=int, required=True, help="Total number of samples in the cohort (e.g., 90).")
    parser.add_argument('--max-af', type=float, default=0.05, help="Maximum allowed internal allele frequency (default: 0.05 / 5%%).")
    parser.add_argument('--annotsv-tsv', help="[Deprecated] External databases (ClinVar/AnnotSV) are not used for AF filtering; AF is calculated strictly internally.")
    parser.add_argument('--output-regions', help="Path to write the unique rare CNV regions (tab-separated CHROM, POS, END) (optional).")
    return parser.parse_args()

def main():
    """
    Main function to calculate variant count and frequency, filter for rare CNVs,
    merge optional AnnotSV annotations, sort the results, and write reports.
    """
    args = parse_args()
    
    if args.total_samples <= 0:
        logging.error(f"Invalid total_samples: {args.total_samples}. Must be a positive integer.")
        sys.exit(1)

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

    # 2. Cluster overlapping CNVs using reciprocal overlap (>= 50%)
    logging.info("Clustering overlapping CNVs across cohort using reciprocal overlap (>= 50%)...")
    df['CLUSTER_ID'] = -1
    cluster_counter = 0
    
    # Sort for deterministic processing
    df_sorted = df.sort_values(by=['CHROM', 'SVTYPE', 'POS'])
    
    groups = df_sorted.groupby(['CHROM', 'SVTYPE'])
    for (chrom, svtype), group_df in groups:
        indices = group_df.index.tolist()
        intervals = group_df[['POS', 'END', 'SAMPLE']].values.tolist()
        
        for i in range(len(intervals)):
            idx_i = indices[i]
            if df.loc[idx_i, 'CLUSTER_ID'] != -1:
                continue
                
            current_cluster = cluster_counter
            cluster_counter += 1
            df.loc[idx_i, 'CLUSTER_ID'] = current_cluster
            
            pos_i, end_i, sample_i = intervals[i]
            len_i = end_i - pos_i
            if len_i <= 0:
                continue
                
            for j in range(i + 1, len(intervals)):
                pos_j, end_j, sample_j = intervals[j]
                if pos_j >= end_i:
                    break
                idx_j = indices[j]
                if df.loc[idx_j, 'CLUSTER_ID'] != -1:
                    continue
                    
                len_j = end_j - pos_j
                if len_j <= 0:
                    continue
                    
                # Check overlap
                overlap_len = min(end_i, end_j) - max(pos_i, pos_j)
                if overlap_len > 0:
                    rec_i = overlap_len / len_i
                    rec_j = overlap_len / len_j
                    # Merge if reciprocal overlap is >= 50%
                    if rec_i >= 0.50 and rec_j >= 0.50:
                        df.loc[idx_j, 'CLUSTER_ID'] = current_cluster

    # 3. Count unique samples per cluster
    cluster_sample_counts = df.groupby('CLUSTER_ID')['SAMPLE'].nunique().reset_index(name='COUNT')
    cluster_sample_counts['INTERNAL_AF'] = cluster_sample_counts['COUNT'] / args.total_samples

    # 4. Map back counts to variants
    df = pd.merge(df, cluster_sample_counts, on='CLUSTER_ID')

    # 5. Filter for Rare CNVs
    final_report = df[df['COUNT'] <= max_allowed_count].copy()
    
    logging.info(f"Total unique clusters before filter: {cluster_counter}")
    logging.info(f"Rare CNV records after filter: {len(final_report)}")

    # Clean report columns
    final_report = final_report.drop(columns=['CLUSTER_ID'])

    # If output-regions is requested, write unique rare CNV coordinates
    if args.output_regions:
        try:
            regions_df = final_report[['CHROM', 'POS', 'END']].drop_duplicates().copy()
            regions_df['CHROM'] = regions_df['CHROM'].astype(str)
            regions_df = regions_df.sort_values(by=['CHROM', 'POS'])
            regions_df.to_csv(args.output_regions, sep='\t', index=False, header=False)
            logging.info(f"Successfully wrote rare CNV regions to {args.output_regions}")
        except Exception as e:
            logging.error(f"Failed to write rare regions file: {e}")

    # Single Responsibility Note: External databases (e.g., ClinVar, gnomAD, AnnotSV)
    # are intentionally NOT merged or consulted here. Allele frequency filtering focuses
    # exclusively on our internal dataset/cohort frequency. External annotations and ACMG
    # classifications are handled separately downstream in addACMGAnnotations.
    if args.annotsv_tsv:
        logging.warning(
            "External annotations (--annotsv-tsv) are ignored during allele frequency filtering. "
            "Allele frequency filtering focuses strictly and exclusively on internal cohort frequency "
            "across our dataset samples, without using external population or clinical databases (e.g. ClinVar, gnomAD). "
            "Downstream ACMG annotations are handled by addACMGAnnotations."
        )

    # 6. Sort for readability (by Chromosome, then Position)
    # Ensure CHROM is treated as string for sorting (e.g., chr1, chr2, chrX)
    final_report['CHROM'] = final_report['CHROM'].astype(str) 
    final_report = final_report.sort_values(by=['CHROM', 'POS'])

    # 7. Save to CSV
    final_report.to_csv(args.output_csv, index=False)
    logging.info(f"Successfully wrote rare CNV report to {args.output_csv}")

if __name__ == '__main__':
    main()
