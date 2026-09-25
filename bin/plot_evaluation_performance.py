#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_evaluation_performance.py
# Purpose: This script generates a multi-panel comparison bar chart evaluating
#   Sensitivity (Recall), Precision (PPV), and Matthews Correlation Coefficient (MCC)
#   between the 'Control' and 'Mosdepth' workflows across different clinical ACMG classes.
# Scope: Cohort evaluation metrics comparison visualization.
# Inputs:
#   Hardcoded data representing performance metrics calculated on validation sets.
# Outputs:
#   /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/evaluate_performance_comparison.png:
#     A 3-panel PNG bar chart showing evaluation metric trends.
# Operational notes:
#   Utilizes a hardcoded data dictionary containing metrics for both workflows
#   under standard classification levels.
# =============================================================================

import os
import matplotlib
matplotlib.use('Agg') # Headless execution
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Set standard plotting style config
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.titlesize': 18
})

def main():
    """
    Constructs the evaluation metrics DataFrame, configures the multi-panel plotting grid,
    and saves the generated comparison bar chart to disk.
    """
    # 1. Prepare clinical evaluation dataset representing performance results
    data = [
        # Control Workflow Evaluation Data
        {"Workflow": "Control", "ACMG": "ACMG 4_5", "Metric": "Sensitivity (%)", "Value": 64.77},
        {"Workflow": "Control", "ACMG": "ACMG 4_5", "Metric": "Precision (%)", "Value": 43.77},
        {"Workflow": "Control", "ACMG": "ACMG 4_5", "Metric": "MCC", "Value": 0.5323},
        
        {"Workflow": "Control", "ACMG": "ACMG 5", "Metric": "Sensitivity (%)", "Value": 64.65},
        {"Workflow": "Control", "ACMG": "ACMG 5", "Metric": "Precision (%)", "Value": 49.23},
        {"Workflow": "Control", "ACMG": "ACMG 5", "Metric": "MCC", "Value": 0.5640},
        
        {"Workflow": "Control", "ACMG": "ACMG 4", "Metric": "Sensitivity (%)", "Value": 0.12},
        {"Workflow": "Control", "ACMG": "ACMG 4", "Metric": "Precision (%)", "Value": 0.74},
        {"Workflow": "Control", "ACMG": "ACMG 4", "Metric": "MCC", "Value": 0.0029},
        
        {"Workflow": "Control", "ACMG": "ACMG 3", "Metric": "Sensitivity (%)", "Value": 16.63},
        {"Workflow": "Control", "ACMG": "ACMG 3", "Metric": "Precision (%)", "Value": 5.88},
        {"Workflow": "Control", "ACMG": "ACMG 3", "Metric": "MCC", "Value": 0.0984},
        
        {"Workflow": "Control", "ACMG": "ACMG 1", "Metric": "Sensitivity (%)", "Value": 1.18},
        {"Workflow": "Control", "ACMG": "ACMG 1", "Metric": "Precision (%)", "Value": 0.84},
        {"Workflow": "Control", "ACMG": "ACMG 1", "Metric": "MCC", "Value": 0.0096},
        
        {"Workflow": "Control", "ACMG": "ACMG NA", "Metric": "Sensitivity (%)", "Value": 0.09},
        {"Workflow": "Control", "ACMG": "ACMG NA", "Metric": "Precision (%)", "Value": 0.05},
        {"Workflow": "Control", "ACMG": "ACMG NA", "Metric": "MCC", "Value": 0.0003},
        
        # Mosdepth Workflow Evaluation Data
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4_5", "Metric": "Sensitivity (%)", "Value": 79.65},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4_5", "Metric": "Precision (%)", "Value": 30.09},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4_5", "Metric": "MCC", "Value": 0.4893},
        
        {"Workflow": "Mosdepth", "ACMG": "ACMG 5", "Metric": "Sensitivity (%)", "Value": 69.61},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 5", "Metric": "Precision (%)", "Value": 28.38},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 5", "Metric": "MCC", "Value": 0.4442},
        
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4", "Metric": "Sensitivity (%)", "Value": 10.04},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4", "Metric": "Precision (%)", "Value": 51.70},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 4", "Metric": "MCC", "Value": 0.2277},
        
        {"Workflow": "Mosdepth", "ACMG": "ACMG 3", "Metric": "Sensitivity (%)", "Value": 2.92},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 3", "Metric": "Precision (%)", "Value": 1.46},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 3", "Metric": "MCC", "Value": 0.0202},
        
        {"Workflow": "Mosdepth", "ACMG": "ACMG 1", "Metric": "Sensitivity (%)", "Value": 0.69},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 1", "Metric": "Precision (%)", "Value": 0.58},
        {"Workflow": "Mosdepth", "ACMG": "ACMG 1", "Metric": "MCC", "Value": 0.0060},
        
        {"Workflow": "Mosdepth", "ACMG": "ACMG NA", "Metric": "Sensitivity (%)", "Value": 0.03},
        {"Workflow": "Mosdepth", "ACMG": "ACMG NA", "Metric": "Precision (%)", "Value": 0.01},
        {"Workflow": "Mosdepth", "ACMG": "ACMG NA", "Metric": "MCC", "Value": -0.0006}
    ]
    
    df = pd.DataFrame(data)
    
    # Create 3 subplots corresponding to: Sensitivity, Precision, and MCC
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=False)
    
    metrics = ["Sensitivity (%)", "Precision (%)", "MCC"]
    titles = [
        "A. Sensitivity (Recall)", 
        "B. Precision (PPV)", 
        "C. Matthews Correlation Coefficient (MCC)"
    ]
    ylabels = ["Sensitivity (%)", "Precision (%)", "Correlation Value"]
    
    # Custom palette colors (Slate Blue for Control vs Dark Coral for Mosdepth)
    colors = {"Control": "#4682B4", "Mosdepth": "#D56C5B"}
    
    for i, metric in enumerate(metrics):
        ax = axes[i]
        sub_df = df[df["Metric"] == metric]
        
        sns.barplot(
            x="ACMG", 
            y="Value", 
            hue="Workflow", 
            data=sub_df, 
            ax=ax, 
            palette=colors,
            edgecolor="black",
            linewidth=1
        )
        
        ax.set_title(titles[i], fontweight="bold", pad=15)
        ax.set_xlabel("ACMG Clinical Prioritization Class", fontweight="bold", labelpad=10)
        ax.set_ylabel(ylabels[i], fontweight="bold", labelpad=10)
        ax.legend(title="Workflow", frameon=True, shadow=False)
        
        # Enable grid overlay
        ax.grid(True, linestyle="--", alpha=0.7)
        
        # Annotate bars with exact decimal values (e.g. 3 decimals for MCC, 1 decimal for %)
        for container in ax.containers:
            fmt = "%.3f" if metric == "MCC" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=4, fontsize=10, weight="semibold")
            
    plt.tight_layout()
    output_path = "/home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/evaluate_performance_comparison.png"
    
    # Ensure parent output directory exists before saving
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    print(f"Successfully generated and saved performance comparison plot to: {output_path}")

if __name__ == "__main__":
    main()
