#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_cnv_qc_metrics.py
# Purpose: This script visualizes call-level quality control (QC) metrics for CNV calling.
#   It loads detailed metrics, identifies attributes like CNV size, GC content, and
#   mappability, and generates a 2x2 diagnostic dashboard panel.
# Scope: Visualization and diagnostic plotting.
# Inputs:
#   --metrics-csv: Path to a CSV file containing call-level QC metrics.
#   --output-dir: Directory path where generated plots will be saved.
# Outputs:
#   cnv_qc_dashboard.png: A 4-panel dashboard containing size distribution,
#     GC content, mappability scores, and calls per sample distributions.
#   no_calls_alert.png: A placeholder image generated if the input call set is empty.
# Operational notes:
#   Ensures robustness by dynamically detecting columns representing CNV type,
#   size, GC content, and mappability.
# =============================================================================

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    """
    Parses command-line arguments.
    
    Returns:
        argparse.Namespace: Object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Visualize call-level CNV QC metrics.")
    parser.add_argument('--metrics-csv', required=True, help="Path to the call-level metrics CSV.")
    parser.add_argument('--output-dir', required=True, help="Directory to save the generated plots.")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. Load call-level metrics dataset
    if not os.path.exists(args.metrics_csv):
        print(f"Error: Metrics file {args.metrics_csv} does not exist.")
        return
        
    try:
        df = pd.read_csv(args.metrics_csv)
    except pd.errors.EmptyDataError:
        # Create an empty DataFrame to trigger the empty dataset handling block below
        df = pd.DataFrame()
    
    # Handle the empty dataset case gracefully by writing a placeholder text plot
    if df.empty:
        print("Warning: Metrics CSV is empty. Generating placeholder alert.")
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No CNVs called to visualize.", ha='center', va='center', fontsize=14)
        plt.savefig(os.path.join(args.output_dir, "cnv_qc_dashboard.png"))
        return

    # Set standard visual theme
    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 14})
    
    # Proactively identify structural variant type column (e.g. DUP vs DEL)
    type_col = 'cnv_type' if 'cnv_type' in df.columns else ('SVTYPE' if 'SVTYPE' in df.columns else None)
    if not type_col:
        # Fallback search strategy: scan column contents for 'DEL' or 'DUP'
        for col in df.columns:
            if df[col].astype(str).str.contains('DEL|DUP').any():
                type_col = col
                break
    
    # Compute size attribute if start and end coordinates are present but size is missing
    if 'size' not in df.columns and 'start' in df.columns and 'end' in df.columns:
        df['size'] = (df['end'] - df['start']).abs()

    # 2. Setup Diagnostic Dashboard Plot (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"CNV Call Quality Control Dashboard ({os.path.basename(args.metrics_csv).split('_')[0]})", fontsize=18, fontweight='bold')

    # Panel A: CNV Size Distribution (Log10 scale to capture small exon-level and large chromosomal events)
    ax_size = axes[0, 0]
    if 'size' in df.columns:
        df['log10_size'] = np.log10(df['size'].replace(0, 1))
        if type_col and df[type_col].nunique() > 1:
            # Color maps for type
            colors = {'DEL': '#e41a1c', 'DUP': '#377eb8'}
            for cat, subdf in df.groupby(type_col):
                subdf['log10_size'].plot(kind='hist', alpha=0.5, ax=ax_size, label=cat, bins=20, color=colors.get(cat, '#4daf4a'))
            ax_size.legend()
        else:
            df['log10_size'].plot(kind='hist', bins=20, ax=ax_size, color="skyblue")
        ax_size.set_title("CNV Size Distribution")
        ax_size.set_xlabel("Log10 Size (bp)")
        ax_size.set_ylabel("Count")
        ax_size.grid(True, linestyle='--', alpha=0.6)
    else:
        ax_size.text(0.5, 0.5, "Size column not found", ha='center', va='center')

    # Panel B: GC Content Distribution of CNV call intervals
    ax_gc = axes[0, 1]
    gc_col = [c for c in df.columns if 'gc' in c.lower()]
    if gc_col:
        gc_name = gc_col[0]
        if type_col and df[type_col].nunique() > 1:
            colors = {'DEL': '#e41a1c', 'DUP': '#377eb8'}
            for cat, subdf in df.groupby(type_col):
                subdf[gc_name].plot(kind='kde', ax=ax_gc, label=cat, color=colors.get(cat, '#4daf4a'))
            ax_gc.legend()
        else:
            df[gc_name].plot(kind='hist', bins=20, ax=ax_gc, color="salmon")
        ax_gc.set_title("GC Content Distribution of CNV Calls")
        ax_gc.set_xlabel("GC Content (%)")
        ax_gc.set_ylabel("Density")
        ax_gc.grid(True, linestyle='--', alpha=0.6)
    else:
        ax_gc.text(0.5, 0.5, "GC Content column not found", ha='center', va='center')

    # Panel C: Mappability Distribution within CNV call intervals (critical for identifying false positives in repeats)
    ax_map = axes[1, 0]
    map_col = [c for c in df.columns if 'map' in c.lower()]
    if map_col:
        map_name = map_col[0]
        if type_col and df[type_col].nunique() > 1:
            colors = {'DEL': '#e41a1c', 'DUP': '#377eb8'}
            for cat, subdf in df.groupby(type_col):
                subdf[map_name].plot(kind='hist', alpha=0.5, ax=ax_map, label=cat, bins=20, color=colors.get(cat, '#4daf4a'))
            ax_map.legend()
        else:
            df[map_name].plot(kind='hist', bins=20, ax=ax_map, color="lightgreen")
        ax_map.set_title("Mappability Score of CNV Calls")
        ax_map.set_xlabel("Mappability")
        ax_map.set_ylabel("Count")
        ax_map.grid(True, linestyle='--', alpha=0.6)
    else:
        ax_map.text(0.5, 0.5, "Mappability column not found", ha='center', va='center')

    # Panel D: Total CNV Call Counts per Sample in the population
    ax_count = axes[1, 1]
    sample_col = 'sample_id' if 'sample_id' in df.columns else ('sample' if 'sample' in df.columns else None)
    if sample_col:
        counts = df[sample_col].value_counts().reset_index()
        counts.columns = ['sample', 'num_cnvs']
        
        # Native Matplotlib boxplot
        ax_count.boxplot(counts['num_cnvs'].values, patch_artist=True, boxprops=dict(facecolor='orchid', alpha=0.6))
        # Native scatter for strip plot (jitter on x-axis)
        x = np.random.normal(1.0, 0.04, size=len(counts))
        ax_count.scatter(x, counts['num_cnvs'].values, color="black", alpha=0.5, zorder=3)
        
        ax_count.set_title("CNV Calls Per Sample")
        ax_count.set_ylabel("Number of Calls")
        ax_count.set_xlabel("Sample Population")
        ax_count.set_xticklabels([])
        ax_count.grid(True, linestyle='--', alpha=0.6)
    else:
        ax_count.text(0.5, 0.5, "Sample ID column not found", ha='center', va='center')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    dashboard_path = os.path.join(args.output_dir, "cnv_qc_dashboard.png")
    plt.savefig(dashboard_path, dpi=150)
    plt.close()
    
    print(f"QC Visualization completed successfully. Dashboard saved to: {dashboard_path}")

if __name__ == '__main__':
    main()
