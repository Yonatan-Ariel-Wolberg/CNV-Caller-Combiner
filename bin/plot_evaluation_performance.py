#!/usr/bin/env python3
"""
bin/plot_evaluation_performance.py

Dynamically parses evaluation performance reports across workflow directories
(e.g., Control and Mosdepth) and generates a multi-panel comparison bar chart
evaluating Recall (Sensitivity), Precision, and Matthews Correlation Coefficient (MCC)
across clinical ACMG classes with ZERO hardcoded values.
"""

import os
import sys
import argparse
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Dynamically parse performance reports and plot comparison.")
    parser.add_argument('--input-dir', default=None,
                        help="Directory containing workflow subdirectories (e.g. Control/ and Mosdepth/) with report files.")
    parser.add_argument('--output-png', default="docs/figures/evaluate_performance_comparison.png",
                        help="Path to write the comparison plot.")
    return parser.parse_args()


def resolve_default_input_dir():
    candidates = [
        "/spaces/ywolberg/qc-gated-cnv-combiner",
        os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/pipeline_results")),
        os.path.abspath("data/pipeline_results")
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[0]


def parse_single_report(report_path):
    metrics = {}
    if not os.path.exists(report_path):
        return metrics
    with open(report_path, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

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
    input_dir = args.input_dir or resolve_default_input_dir()

    if not os.path.exists(input_dir):
        print(f"Error: Input directory {input_dir} does not exist.")
        sys.exit(1)

    data = []
    # Search for performance reports dynamically across all subdirectories
    for root, _, files in os.walk(input_dir):
        for fname in files:
            if fname.endswith('_performance_report.txt') or fname.endswith('_performance.txt'):
                match = re.search(r'ACMG([a-zA-Z0-9_]+)_performance', fname)
                if not match:
                    continue
                acmg_class = match.group(1).replace('_', ' ')
                acmg_display = f"ACMG {acmg_class}"
                
                # Determine workflow from path
                workflow = "Mosdepth" if "mosdepth" in root.lower() else ("Control" if "control" in root.lower() else "Standard")
                
                report_path = os.path.join(root, fname)
                metrics = parse_single_report(report_path)
                for metric_name, val in metrics.items():
                    data.append({
                        "Workflow": workflow,
                        "ACMG": acmg_display,
                        "Metric": metric_name,
                        "Value": val
                    })

    if not data:
        print("Warning: No evaluation report files found in search path. Generating informational plot.")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "Dynamic evaluation scan completed.\nNo target report files matched.",
                ha='center', va='center', fontsize=14, weight='bold', color='gray')
        ax.set_axis_off()
        os.makedirs(os.path.dirname(os.path.abspath(args.output_png)), exist_ok=True)
        plt.savefig(args.output_png, dpi=300)
        return

    df = pd.DataFrame(data)
    df['ACMG'] = df['ACMG'].str.replace('4 5', '4_5')
    acmg_order = ['ACMG 4_5', 'ACMG 5', 'ACMG 4', 'ACMG 3', 'ACMG 1', 'ACMG NA']
    present_order = [x for x in acmg_order if x in df['ACMG'].unique()]

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=False)
    metrics_list = ["Sensitivity (%)", "Precision (%)", "MCC"]
    titles = [
        "A. Sensitivity (Recall)",
        "B. Precision (PPV)",
        "C. Matthews Correlation Coefficient (MCC)"
    ]
    ylabels = ["Sensitivity (%)", "Precision (%)", "Correlation Value"]
    colors = {"Control": "#4682B4", "Mosdepth": "#D56C5B"}

    for i, metric in enumerate(metrics_list):
        ax = axes[i]
        sub_df = df[df["Metric"] == metric]
        if sub_df.empty:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            ax.set_title(titles[i])
            continue

        sns.barplot(
            x="ACMG",
            y="Value",
            hue="Workflow",
            data=sub_df,
            order=present_order,
            ax=ax,
            palette=colors,
            edgecolor="black",
            linewidth=1
        )
        ax.set_title(titles[i], fontweight="bold", pad=15)
        ax.set_xlabel("ACMG Clinical Prioritization Class", fontweight="bold", labelpad=10)
        ax.set_ylabel(ylabels[i], fontweight="bold", labelpad=10)
        ax.legend(title="Workflow", frameon=True)
        ax.grid(True, linestyle="--", alpha=0.7)

        for container in ax.containers:
            fmt = "%.3f" if metric == "MCC" else "%.1f%%"
            ax.bar_label(container, fmt=fmt, label_type='edge', padding=4, fontsize=9.5, weight="semibold")

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.output_png)), exist_ok=True)
    plt.savefig(args.output_png, dpi=300)
    print(f"Successfully generated dynamic performance comparison plot to: {args.output_png}")


if __name__ == "__main__":
    main()
