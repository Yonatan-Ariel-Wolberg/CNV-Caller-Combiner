#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_combination_performance.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Headless run
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set standard plotting style
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.titlesize': 16
})

def main():
    csv_path = "results_e2e_comparison/caller_combinations_acmg_metrics.csv"
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)
    
    # Let's focus on the key combinations and categories
    key_combos = [
        "dragen + canoes",
        "dragen + xhmm",
        "dragen + gatk-gcnv",
        "canoes + xhmm",
        "cnvkit + canoes"
    ]
    
    df_filtered = df[df["Combination"].isin(key_combos)].copy()
    
    # Convert sensitivity and precision to percentages for plotting
    df_filtered["Sensitivity_Pct"] = df_filtered["Sensitivity"] * 100
    df_filtered["Precision_Pct"] = df_filtered["Precision"] * 100
    
    # ----------------------------------------------------
    # Plot 1: ACMG 5 (Pathogenic) Performance Comparison
    # ----------------------------------------------------
    df_acmg5 = df_filtered[df_filtered["Category"] == "ACMG5"].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    metrics = ["Sensitivity_Pct", "Precision_Pct", "F_beta"]
    ylabels = ["Sensitivity (%)", "Precision (%)", "F2-Score"]
    titles = ["A. Sensitivity (Recall) - ACMG5", "B. Precision (PPV) - ACMG5", "C. F2-Score - ACMG5"]
    
    colors = {"Union": "#D56C5B", "Intersection": "#4682B4"}
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        sns.barplot(
            x="Combination",
            y=metric,
            hue="Type",
            data=df_acmg5,
            ax=ax,
            palette=colors,
            edgecolor="black",
            linewidth=1
        )
        ax.set_title(titles[idx], fontweight="bold", pad=12)
        ax.set_ylabel(ylabels[idx], fontweight="bold")
        ax.set_xlabel("Caller Combination", fontweight="bold")
        ax.tick_params(axis='x', rotation=15)
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(title="Combination Logic")
        
        # Annotate bars
        for container in ax.containers:
            fmt = "%.3f" if metric == "F_beta" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=3, fontsize=9, weight="semibold")
            
    plt.tight_layout()
    output_acmg5 = "results_e2e_comparison/caller_combinations_acmg5_performance.png"
    plt.savefig(output_acmg5, dpi=300)
    print(f"Generated plot: {output_acmg5}")
    
    # ----------------------------------------------------
    # Plot 2: Total (Overall Cohort) Performance Comparison
    # ----------------------------------------------------
    df_total = df_filtered[df_filtered["Category"] == "Total"].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    titles_total = ["A. Sensitivity (Recall) - Overall", "B. Precision (PPV) - Overall", "C. F2-Score - Overall"]
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        sns.barplot(
            x="Combination",
            y=metric,
            hue="Type",
            data=df_total,
            ax=ax,
            palette=colors,
            edgecolor="black",
            linewidth=1
        )
        ax.set_title(titles_total[idx], fontweight="bold", pad=12)
        ax.set_ylabel(ylabels[idx], fontweight="bold")
        ax.set_xlabel("Caller Combination", fontweight="bold")
        ax.tick_params(axis='x', rotation=15)
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(title="Combination Logic")
        
        # Annotate bars
        for container in ax.containers:
            fmt = "%.3f" if metric == "F_beta" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=3, fontsize=9, weight="semibold")
            
    plt.tight_layout()
    output_total = "results_e2e_comparison/caller_combinations_total_performance.png"
    plt.savefig(output_total, dpi=300)
    print(f"Generated plot: {output_total}")

if __name__ == "__main__":
    main()
