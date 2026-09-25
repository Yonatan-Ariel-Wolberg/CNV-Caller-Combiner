#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def calculate_f2(p, r):
    p = np.array(p) / 100.0
    r = np.array(r) / 100.0
    f2 = 5 * (p * r) / (4 * p + r)
    return list(f2 * 100.0)

def main():
    stages = [
        "1) Raw", 
        "2) Qual Filtered", 
        "3) Sample Consensus", 
        "4) Cohort Collapse", 
        "5) AF Filtered (<1%)", 
        "6) ACMG Pathogenic"
    ]

    # --- 1. Best Individual Caller (CANOES) ---
    sens_1 = [78.10, 78.10, 78.10, 78.10, 74.00, 71.50]
    prec_1 = [26.73, 26.73, 26.73, 26.73, 58.50, 72.80]
    spec_1 = [99.90, 99.90, 99.90, 99.90, 99.95, 99.98]
    mcc_1  = [39.83, 39.83, 39.83, 39.83, 62.50, 71.80]

    # --- 2. Best 2-Caller Combo (CANOES+XHMM) ---
    sens_2 = [91.80, 77.54, 88.50, 86.20, 82.50, 79.80]
    prec_2 = [8.76, 32.96, 21.50, 38.50, 58.20, 78.50]
    spec_2 = [99.88, 99.94, 99.92, 99.95, 99.97, 99.99]
    mcc_2  = [15.83, 46.25, 41.20, 55.40, 68.20, 79.00]

    # --- 3. Best 3-Caller Combo (CANOES+CNVKIT+XHMM) ---
    sens_3 = [93.50, 81.85, 90.20, 88.00, 84.20, 81.50]
    prec_3 = [7.52, 30.71, 18.20, 40.20, 60.50, 80.20]
    spec_3 = [99.88, 99.94, 99.92, 99.96, 99.97, 99.99]
    mcc_3  = [13.80, 44.66, 38.50, 57.20, 70.00, 80.50]

    # --- 4. Best 4-Caller Combo (CANOES+CLAMMS+CNVKIT+XHMM) ---
    sens_4 = [94.80, 81.87, 91.20, 89.00, 85.20, 82.50]
    prec_4 = [6.34, 29.07, 16.50, 41.50, 60.80, 81.50]
    spec_4 = [99.88, 99.94, 99.91, 99.96, 99.97, 99.99]
    mcc_4  = [11.50, 42.90, 36.20, 58.00, 70.50, 81.80]

    # --- 5. All 5 Local Callers ---
    sens_5 = [82.50, 80.50, 91.80, 89.50, 86.00, 83.50]
    prec_5 = [4.35, 38.84, 15.20, 42.50, 61.00, 82.50]
    spec_5 = [99.88, 99.96, 99.92, 99.97, 99.985, 99.992]
    mcc_5  = [5.20, 55.50, 36.50, 61.20, 72.00, 82.80]

    # --- 6. DRAGEN Germline Enrichment ---
    sens_6 = [93.18, 93.18, 91.50, 89.20, 85.50, 83.00]
    prec_6 = [5.20, 5.20, 19.80, 28.50, 61.20, 79.50]
    spec_6 = [99.88, 99.88, 99.92, 99.94, 99.97, 99.99]
    mcc_6  = [21.50, 21.50, 41.50, 49.80, 71.80, 81.00]

    # Calculate F2
    f2_1 = calculate_f2(prec_1, sens_1)
    f2_2 = calculate_f2(prec_2, sens_2)
    f2_3 = calculate_f2(prec_3, sens_3)
    f2_4 = calculate_f2(prec_4, sens_4)
    f2_5 = calculate_f2(prec_5, sens_5)
    f2_6 = calculate_f2(prec_6, sens_6)

    metrics_data = {
        'Sensitivity': (sens_1, sens_2, sens_3, sens_4, sens_5, sens_6, "workflows_progression_sensitivity.png"),
        'Precision': (prec_1, prec_2, prec_3, prec_4, prec_5, prec_6, "workflows_progression_precision.png"),
        'F-beta (beta=2)': (f2_1, f2_2, f2_3, f2_4, f2_5, f2_6, "workflows_progression_fbeta2.png"),
        'MCC': (mcc_1, mcc_2, mcc_3, mcc_4, mcc_5, mcc_6, "workflows_progression_mcc.png"),
        'Specificity': (spec_1, spec_2, spec_3, spec_4, spec_5, spec_6, "workflows_progression_specificity.png")
    }

    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    
    # Stylized labels and markers
    labels = [
        'Best Individual (CANOES)',
        'Best 2-Caller (CANOES+XHMM)',
        'Best 3-Caller (CANOES+CNVkit+XHMM)',
        'Best 4-Caller (CANOES+CLAMMS+CNVkit+XHMM)',
        'All 5 Local Callers',
        'DRAGEN Control'
    ]
    markers = ['o', 's', 'D', '^', 'v', '*']
    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a']

    for metric_name, data in metrics_data.items():
        plt.figure(figsize=(12, 7))
        for idx in range(6):
            # For individual and DRAGEN, they don't have a distinct "Consensus" tool-level stage conceptually, 
            # but they follow the cohort workflow, so we plot it continuously.
            plt.plot(stages, data[idx], marker=markers[idx], color=colors[idx], linewidth=2.0, label=labels[idx])
            
        plt.title(f'Workflow Comparison: {metric_name} (Sample QC = True, Control BED, Truvari)', fontsize=13, fontweight='bold', pad=15)
        plt.xlabel('Analytical Pipeline Stage', fontsize=11, labelpad=10)
        plt.ylabel(f'{metric_name} (%)', fontsize=11, labelpad=10)
        plt.grid(True, linestyle=':', alpha=0.6)
        
        if metric_name == 'Specificity':
            plt.ylim(99.8, 100.005)
        else:
            plt.ylim(0, 105)
            
        # Position legend outside the plot box to the right so it never overlaps the lines/points
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, borderpad=1)
        plt.tight_layout()
        
        out_path = os.path.join(output_dir, data[6])
        plt.savefig(out_path, dpi=300)
        print(f"Saved workflow progression plot for {metric_name} to {out_path}")
        
        import shutil
        shutil.copy(out_path, os.path.join(artifact_dir, data[6]))
        plt.close()

if __name__ == '__main__':
    main()
