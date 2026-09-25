#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: plot_evaluation_performance_dynamic.py
# Purpose: This script scans a target directory containing workflow-specific folders
#   (e.g., Control/ and Mosdepth/), dynamically parses text evaluation report files
#   using regular expressions, compiles metrics across ACMG priority groups, and
#   generates a comparative 3-panel bar chart.
# Scope: Dynamic evaluation report parsing and comparative visualization.
# Inputs:
#   --input-dir: Directory containing subdirectories (workflow runs) with text reports.
#   --output-png: Output image file path for the comparative bar plot.
# Outputs:
#   --output-png: A 3-panel PNG plot (Sensitivity, Precision, and MCC).
# Operational notes:
#   Gracefully handles missing directories or reports by producing a warning placeholder image.
# =============================================================================

import os
import argparse
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Headless run
import matplotlib.pyplot as plt

def parse_args():
    """
    Parses command-line arguments.
    
    Returns:
        argparse.Namespace: Object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Parse performance reports dynamically and plot comparison.")
    parser.add_argument('--input-dir', required=True, help="Directory containing subdirectories (e.g. Control/ and Mosdepth/) with report files.")
    parser.add_argument('--output-png', required=True, help="Path to write the comparison plot.")
    return parser.parse_args()

def parse_report(report_path):
    """
    Parses a single performance report text file using regex to extract Sensitivity,
    Precision, and Matthews Correlation Coefficient (MCC).
    
    Args:
        report_path (str): Absolute or relative path to the report text file.
        
    Returns:
        dict: Extracted metrics keys paired with float values.
    """
    metrics = {}
    with open(report_path, 'r') as f:
        content = f.read()
        
        # Regular expressions to locate metric patterns in text report lines
        sens_match = re.search(r'Sensitivity \(Recall\)\s*:\s*([0-9\.]+)', content)
        prec_match = re.search(r'Precision \(PPV\)\s*:\s*([0-9\.]+)', content)
        mcc_match = re.search(r'Matthews Corr Coeff\s*:\s*([\-\+0-9\.]+)', content)
        
        if sens_match:
            metrics['Sensitivity (%)'] = float(sens_match.group(1)) * 100
        if prec_match:
            metrics['Precision (%)'] = float(prec_match.group(1)) * 100
        if mcc_match:
            metrics['MCC'] = float(mcc_match.group(1))
            
    return metrics

def main():
    args = parse_args()
    data = []
    
    # 1. Validation of the base report directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist.")
        return
        
    # Walk workflow subfolders (e.g. Control, Mosdepth)
    for sub in os.listdir(args.input_dir):
        sub_path = os.path.join(args.input_dir, sub)
        if not os.path.isdir(sub_path):
            continue
            
        workflow = sub
        # Find ACMG performance reports in folder
        for fname in os.listdir(sub_path):
            if fname.endswith('_performance_report.txt'):
                match = re.search(r'ACMG([a-zA-Z0-9_]+)_performance_report\.txt', fname)
                if not match:
                    continue
                # Normalize clinical classification tags (e.g. 4_5)
                acmg_class = match.group(1).replace('_', ' ')
                acmg_display = f"ACMG {acmg_class}"
                
                report_path = os.path.join(sub_path, fname)
                metrics = parse_report(report_path)
                
                for metric_name, val in metrics.items():
                    data.append({
                        "Workflow": workflow,
                        "ACMG": acmg_display,
                        "Metric": metric_name,
                        "Value": val
                    })
                    
    # Handle empty/missing data scenario by creating an elegant placeholder plot
    if not data:
        print("No evaluation report files found or successfully parsed. Generating blank comparison plot.")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No variants called or evaluated\nin the current run context.", 
                ha='center', va='center', fontsize=14, weight='bold', color='gray')
        ax.set_axis_off()
        os.makedirs(os.path.dirname(os.path.abspath(args.output_png)), exist_ok=True)
        plt.savefig(args.output_png, dpi=300)
        return
        
    df = pd.DataFrame(data)
    df['ACMG'] = df['ACMG'].str.replace('4 5', '4_5')
    
    # Enforce priority categorical ordering for clinical classifications
    acmg_order = ['ACMG 4_5', 'ACMG 5', 'ACMG 4', 'ACMG 3', 'ACMG 1', 'ACMG NA']
    present_order = [x for x in acmg_order if x in df['ACMG'].unique()]
    
    # Configure plotting parameters
    plt.style.use('ggplot')
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.titlesize': 18
    })
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=False)
    metrics_list = ["Sensitivity (%)", "Precision (%)", "MCC"]
    titles = [
        "A. Sensitivity (Recall)", 
        "B. Precision (PPV)", 
        "C. Matthews Correlation Coefficient (MCC)"
    ]
    ylabels = ["Sensitivity (%)", "Precision (%)", "Correlation Value"]
    
    # Establish a color scheme matching default workflows or falling back dynamically
    workflows_present = df['Workflow'].unique()
    colors_default = ["#4682B4", "#D56C5B", "#2CA02C", "#9467BD", "#8C564B"]
    colors = {}
    for idx, wf in enumerate(workflows_present):
        if wf == "Control":
            colors[wf] = "#4682B4"
        elif wf == "Mosdepth":
            colors[wf] = "#D56C5B"
        else:
            colors[wf] = colors_default[idx % len(colors_default)]
            
    # Draw bars across panels
    for i, metric in enumerate(metrics_list):
        ax = axes[i]
        sub_df = df[df["Metric"] == metric]
        
        if sub_df.empty:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            ax.set_title(titles[i])
            continue
            
        pivot_df = sub_df.pivot(index='ACMG', columns='Workflow', values='Value').reindex(present_order)
        pivot_df.plot(
            kind='bar',
            ax=ax,
            color=[colors[col] for col in pivot_df.columns if col in colors],
            edgecolor="black",
            linewidth=1,
            rot=0
        )
        
        ax.set_title(titles[i], fontweight="bold", pad=15)
        ax.set_xlabel("ACMG Clinical Prioritization Class", fontweight="bold", labelpad=10)
        ax.set_ylabel(ylabels[i], fontweight="bold", labelpad=10)
        ax.legend(title="Workflow", frameon=True, shadow=False)
        ax.grid(True, linestyle="--", alpha=0.7)
        
        # Label bars with exact values
        for container in ax.containers:
            fmt = "%.3f" if metric == "MCC" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=4, fontsize=10, weight="semibold")
            
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.output_png)), exist_ok=True)
    plt.savefig(args.output_png, dpi=300)
    print(f"Successfully generated and saved comparison plot to: {args.output_png}")

if __name__ == '__main__':
    main()
