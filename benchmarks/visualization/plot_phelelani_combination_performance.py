#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_phelelani_combination_performance.py
# Purpose: This script visualizes caller combinations performance metrics (Sensitivity,
#   Precision, and F2-Score) on a validation dataset, evaluating Union and Intersection
#   combination styles for combinations published by Phelelani.
# Scope: Visualization of published pipeline validation metrics.
# Inputs:
#   results_e2e_comparison/phelelani_caller_combinations_metrics.csv: CSV listing metrics
#     representing Sensitivity, Precision, and F-beta scores per caller combination.
# Outputs:
#   results_e2e_comparison/phelelani_caller_combinations_performance.png:
#     A 3-panel comparison bar chart detailing performance characteristics.
# Operational notes:
#   Supports quick validation comparisons.
# =============================================================================

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Headless run
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set standard plotting style config
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
    """
    Loads combination metrics from the validation run, formats metric columns,
    and plots bar charts comparing Union and Intersection performance metrics.
    """
    csv_path = "results_e2e_comparison/phelelani_caller_combinations_metrics.csv"
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)
    
    # Convert sensitivity and precision to percentages for plotting
    df["Sensitivity_Pct"] = df["Sensitivity"] * 100
    df["Precision_Pct"] = df["Precision"] * 100
    
    # Create 3 subplots corresponding to: Sensitivity, Precision, and F2-Score
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    metrics = ["Sensitivity_Pct", "Precision_Pct", "F_beta"]
    ylabels = ["Sensitivity (%)", "Precision (%)", "F2-Score"]
    titles = ["A. Sensitivity (Recall)", "B. Precision (PPV)", "C. F2-Score"]
    
    # Establish a color scheme matching Union vs Intersection logics
    colors = {"Union": "#D56C5B", "Intersection": "#4682B4"}
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        sns.barplot(
            x="Combination",
            y=metric,
            hue="Type",
            data=df,
            ax=ax,
            palette=colors,
            edgecolor="black",
            linewidth=1
        )
        ax.set_title(titles[idx], fontweight="bold", pad=12)
        ax.set_ylabel(ylabels[idx], fontweight="bold")
        ax.set_xlabel("Caller Combination", fontweight="bold")
        ax.tick_params(axis='x', rotation=10)
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(title="Combination Logic")
        
        # Annotate bars with numeric labels
        for container in ax.containers:
            fmt = "%.3f" if metric == "F_beta" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=3, fontsize=9, weight="semibold")
            
    plt.tight_layout()
    output_png = "results_e2e_comparison/phelelani_caller_combinations_performance.png"
    plt.savefig(output_png, dpi=300)
    print(f"Generated plot: {output_png}")

if __name__ == "__main__":
    main()
