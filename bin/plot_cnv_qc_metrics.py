#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_cnv_qc_metrics.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize call-level CNV QC metrics.")
    parser.add_argument('--metrics-csv', required=True, help="Path to the call-level metrics CSV.")
    parser.add_argument('--output-dir', required=True, help="Directory to save the generated plots.")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. Load Data
    if not os.path.exists(args.metrics_csv):
        print(f"Error: Metrics file {args.metrics_csv} does not exist.")
        return
        
    df = pd.read_csv(args.metrics_csv)
    if df.empty:
        print("Warning: Metrics CSV is empty. Generating placeholder alert.")
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No CNVs called to visualize.", ha='center', va='center', fontsize=14)
        plt.savefig(os.path.join(args.output_dir, "no_calls_alert.png"))
        return

    # Set styling
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 14})
    
    # Standardize CNV Type Column name if needed
    type_col = 'cnv_type' if 'cnv_type' in df.columns else ('SVTYPE' if 'SVTYPE' in df.columns else None)
    if not type_col:
        # Fallback to identify column with DEL/DUP
        for col in df.columns:
            if df[col].astype(str).str.contains('DEL|DUP').any():
                type_col = col
                break
    
    # Add size column if start/end exist
    if 'size' not in df.columns and 'start' in df.columns and 'end' in df.columns:
        df['size'] = (df['end'] - df['start']).abs()

    # 2. Setup Dashboard Plot (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"CNV Call Quality Control Dashboard ({os.path.basename(args.metrics_csv).split('_')[0]})", fontsize=18, fontweight='bold')

    # A. Size Distribution
    ax_size = axes[0, 0]
    if 'size' in df.columns:
        # log10 scale for CNV sizes
        df['log10_size'] = np.log10(df['size'].replace(0, 1))
        if type_col and df[type_col].nunique() > 1:
            sns.histplot(data=df, x='log10_size', hue=type_col, kde=True, ax=ax_size, multiple="dodge", palette="Set1")
        else:
            sns.histplot(data=df, x='log10_size', kde=True, ax=ax_size, color="skyblue")
        ax_size.set_title("CNV Size Distribution")
        ax_size.set_xlabel("Log10 Size (bp)")
        ax_size.set_ylabel("Count")
    else:
        ax_size.text(0.5, 0.5, "Size column not found", ha='center', va='center')

    # B. GC Content Distribution
    ax_gc = axes[0, 1]
    gc_col = [c for c in df.columns if 'gc' in c.lower()]
    if gc_col:
        gc_name = gc_col[0]
        if type_col and df[type_col].nunique() > 1:
            sns.kdeplot(data=df, x=gc_name, hue=type_col, fill=True, ax=ax_gc, common_norm=False, palette="Set1")
        else:
            sns.histplot(data=df, x=gc_name, kde=True, ax=ax_gc, color="salmon")
        ax_gc.set_title("GC Content Distribution of CNV Calls")
        ax_gc.set_xlabel("GC Content (%)")
        ax_gc.set_ylabel("Density")
    else:
        ax_gc.text(0.5, 0.5, "GC Content column not found", ha='center', va='center')

    # C. Mappability Distribution
    ax_map = axes[1, 0]
    map_col = [c for c in df.columns if 'map' in c.lower()]
    if map_col:
        map_name = map_col[0]
        if type_col and df[type_col].nunique() > 1:
            sns.histplot(data=df, x=map_name, hue=type_col, multiple="stack", ax=ax_map, palette="Set1")
        else:
            sns.histplot(data=df, x=map_name, kde=True, ax=ax_map, color="lightgreen")
        ax_map.set_title("Mappability Score of CNV Calls")
        ax_map.set_xlabel("Mappability")
        ax_map.set_ylabel("Count")
    else:
        ax_map.text(0.5, 0.5, "Mappability column not found", ha='center', va='center')

    # D. CNV Count Per Sample
    ax_count = axes[1, 1]
    sample_col = 'sample_id' if 'sample_id' in df.columns else ('sample' if 'sample' in df.columns else None)
    if sample_col:
        counts = df[sample_col].value_counts().reset_index()
        counts.columns = ['sample', 'num_cnvs']
        sns.boxplot(data=counts, y='num_cnvs', ax=ax_count, color="orchid")
        sns.stripplot(data=counts, y='num_cnvs', ax=ax_count, color="black", alpha=0.5, jitter=0.2)
        ax_count.set_title("CNV Calls Per Sample")
        ax_count.set_ylabel("Number of Calls")
        ax_count.set_xlabel("Sample Population")
    else:
        ax_count.text(0.5, 0.5, "Sample ID column not found", ha='center', va='center')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    dashboard_path = os.path.join(args.output_dir, "cnv_qc_dashboard.png")
    plt.savefig(dashboard_path, dpi=150)
    plt.close()
    
    print(f"QC Visualization completed successfully. Dashboard saved to: {dashboard_path}")

if __name__ == '__main__':
    main()
