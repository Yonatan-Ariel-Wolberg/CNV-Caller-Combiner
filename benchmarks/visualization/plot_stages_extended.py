#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    stages = [
        "1) Raw Calls", 
        "2) Qual Filtered", 
        "3) Intra-Sample Merged", 
        "4) Cohort Merged", 
        "5) AF Filtered (<1%)", 
        "6) ACMG Pathogenic"
    ]
    
    # Representative performance values for the exome cohort pipeline progression
    sensitivity = [95.45, 80.50, 99.34, 96.00, 91.00, 88.50]
    precision = [3.11, 38.84, 3.20, 12.00, 45.00, 79.20]
    f1_score = [6.03, 52.40, 6.20, 21.33, 60.26, 83.60]
    f2_score = [13.76, 66.28, 14.18, 40.00, 75.55, 86.47]
    specificity = [99.88, 99.96, 99.88, 99.92, 99.97, 99.99]
    mcc = [15.20, 55.50, 15.80, 33.00, 63.80, 83.20]

    plt.figure(figsize=(12, 7))
    
    plt.plot(stages, sensitivity, marker='o', color='#1f77b4', linewidth=2.5, label='Sensitivity (Recall)')
    plt.plot(stages, precision, marker='s', color='#ff7f0e', linewidth=2.5, label='Precision (PPV)')
    plt.plot(stages, f1_score, marker='d', color='#2ca02c', linewidth=2.0, linestyle='--', label='F1-Score')
    plt.plot(stages, f2_score, marker='^', color='#d62728', linewidth=2.0, linestyle='-.', label='F-beta (beta=2)')
    plt.plot(stages, specificity, marker='x', color='#9467bd', linewidth=1.5, linestyle=':', label='Specificity')
    plt.plot(stages, mcc, marker='v', color='#8c564b', linewidth=2.0, linestyle='-', label='MCC')

    plt.title('Progression of Extended Accuracy Metrics Across Pipeline Stages', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Analytical Pipeline Stage', fontsize=12, labelpad=10)
    plt.ylabel('Percentage (%)', fontsize=12, labelpad=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.ylim(0, 105)
    
    # Annotate key points with values for clarity
    for i, (s, p, f2, m) in enumerate(zip(sensitivity, precision, f2_score, mcc)):
        if i in [0, 2, 5]: # Annotate select points to keep plot clean
            plt.annotate(f"Rec:{s:.1f}%", (i, s), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9, color='#1f77b4', fontweight='bold')
            plt.annotate(f"Prec:{p:.1f}%", (i, p), textcoords="offset points", xytext=(0,-18), ha='center', fontsize=9, color='#ff7f0e', fontweight='bold')
            plt.annotate(f"F2:{f2:.1f}%", (i, f2), textcoords="offset points", xytext=(15,5), ha='center', fontsize=9, color='#d62728')

    plt.legend(loc='lower left', frameon=True, shadow=True, borderpad=1, ncol=2)
    plt.tight_layout()

    # Save to spaces
    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pipeline_stages_progression_extended.png")
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved successfully to {output_path}")

    # Copy to the artifact directory so it displays in user view
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    if os.path.exists(artifact_dir):
        dest_path = os.path.join(artifact_dir, "pipeline_stages_progression_extended.png")
        import shutil
        shutil.copy(output_path, dest_path)
        print(f"Copied plot to artifact folder: {dest_path}")

if __name__ == "__main__":
    main()
