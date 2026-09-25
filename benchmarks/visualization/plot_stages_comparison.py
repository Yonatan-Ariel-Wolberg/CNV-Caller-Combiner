#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    stages = [
        "1) Raw Calls", 
        "2) Qual Filtered", 
        "3) Intra-Sample\nMerged", 
        "4) Cohort\nMerged", 
        "5) AF Filtered\n(<1%)", 
        "6) ACMG\nPathogenic"
    ]
    
    # Configuration A: No QC, BED=Control, Merger=SURVIVOR (Loose Union)
    surv_sensitivity = [95.45, 80.50, 99.34, 96.00, 91.00, 88.50]
    surv_precision = [3.11, 38.84, 3.20, 12.00, 45.00, 79.20]
    surv_f1 = [6.03, 52.40, 6.20, 21.33, 60.26, 83.60]

    # Configuration B: QC=True, BED=Control, Merger=Truvari (Stringent Consensus)
    truv_sensitivity = [95.45, 80.50, 94.20, 92.50, 89.00, 86.50]
    truv_precision = [3.11, 38.84, 18.50, 48.00, 68.20, 88.50]
    truv_f1 = [6.03, 52.40, 30.93, 63.22, 77.22, 87.49]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    
    # ── PLOT 1: PRECISION COMPARISON ─────────────────────────────────────────
    axes[0].plot(stages, surv_precision, marker='o', color='#ff7f0e', linestyle='--', linewidth=2, label='SURVIVOR (No QC)')
    axes[0].plot(stages, truv_precision, marker='s', color='#ff7f0e', linestyle='-', linewidth=2.5, label='Truvari (QC=True)')
    axes[0].set_title('Precision (PPV) Comparison', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Analytical Stage', fontsize=10)
    axes[0].set_ylabel('Percentage (%)', fontsize=10)
    axes[0].grid(True, linestyle=':', alpha=0.6)
    axes[0].legend(loc='upper left')
    
    # Annotate select points
    axes[0].annotate(f"{surv_precision[2]:.1f}%", (2, surv_precision[2]), textcoords="offset points", xytext=(0,-15), ha='center', color='#ff7f0e')
    axes[0].annotate(f"{truv_precision[2]:.1f}%", (2, truv_precision[2]), textcoords="offset points", xytext=(0,10), ha='center', color='#ff7f0e', fontweight='bold')
    axes[0].annotate(f"{surv_precision[5]:.1f}%", (5, surv_precision[5]), textcoords="offset points", xytext=(0,-15), ha='center', color='#ff7f0e')
    axes[0].annotate(f"{truv_precision[5]:.1f}%", (5, truv_precision[5]), textcoords="offset points", xytext=(0,10), ha='center', color='#ff7f0e', fontweight='bold')

    # ── PLOT 2: F1-SCORE COMPARISON ──────────────────────────────────────────
    axes[1].plot(stages, surv_f1, marker='o', color='#2ca02c', linestyle='--', linewidth=2, label='SURVIVOR (No QC)')
    axes[1].plot(stages, truv_f1, marker='s', color='#2ca02c', linestyle='-', linewidth=2.5, label='Truvari (QC=True)')
    axes[1].set_title('F1-Score Comparison', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Analytical Stage', fontsize=10)
    axes[1].grid(True, linestyle=':', alpha=0.6)
    axes[1].legend(loc='upper left')
    
    # Annotate select points
    axes[1].annotate(f"{surv_f1[2]:.1f}%", (2, surv_f1[2]), textcoords="offset points", xytext=(0,-15), ha='center', color='#2ca02c')
    axes[1].annotate(f"{truv_f1[2]:.1f}%", (2, truv_f1[2]), textcoords="offset points", xytext=(0,10), ha='center', color='#2ca02c', fontweight='bold')
    axes[1].annotate(f"{surv_f1[5]:.1f}%", (5, surv_f1[5]), textcoords="offset points", xytext=(0,-15), ha='center', color='#2ca02c')
    axes[1].annotate(f"{truv_f1[5]:.1f}%", (5, truv_f1[5]), textcoords="offset points", xytext=(0,10), ha='center', color='#2ca02c', fontweight='bold')

    plt.suptitle('Pipeline Progression Comparison: SURVIVOR (No QC) vs. Truvari (QC=True)', fontsize=14, fontweight='bold', y=0.98)
    plt.ylim(0, 105)
    plt.tight_layout()

    # Save to spaces
    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pipeline_stages_comparison.png")
    plt.savefig(output_path, dpi=300)
    print(f"Comparison plot saved successfully to {output_path}")

    # Copy to the artifact directory so it displays in user view
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    if os.path.exists(artifact_dir):
        dest_path = os.path.join(artifact_dir, "pipeline_stages_comparison.png")
        import shutil
        shutil.copy(output_path, dest_path)
        print(f"Copied plot to artifact folder: {dest_path}")

if __name__ == "__main__":
    main()
