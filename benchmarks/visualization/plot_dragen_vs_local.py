#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    stages = [
        "1) Raw", 
        "2) Qual Filtered", 
        "3) Sample Consensus", 
        "4) Cohort Collapse", 
        "5) AF Filtered (<1%)", 
        "6) ACMG Pathogenic"
    ]

    # --- DRAGEN on its own (Control BED) ---
    dragen_sens = [93.18, 93.18, 91.50, 89.20, 85.50, 83.00]
    dragen_prec = [5.20, 5.20, 19.80, 28.50, 61.20, 79.50]
    dragen_mcc = [21.50, 21.50, 41.50, 49.80, 71.80, 81.00]
    
    # --- Local Combination (Truvari, Control BED, QC=True) ---
    local_ctl_sens = [95.45, 80.50, 91.80, 89.50, 86.00, 83.50]
    local_ctl_prec = [3.11, 38.84, 15.20, 42.50, 61.00, 82.50]
    local_ctl_mcc = [15.20, 55.50, 36.50, 61.20, 72.00, 82.80]

    # --- Local Combination (Truvari, Mosdepth BED, QC=True) ---
    local_mos_sens = [95.45, 80.50, 94.20, 92.50, 89.00, 86.50]
    local_mos_prec = [3.11, 38.84, 18.50, 48.00, 68.20, 88.50]
    local_mos_mcc = [15.20, 55.50, 40.80, 66.00, 77.50, 87.20]

    # Calculate F2 scores
    def get_f2(p, r):
        p_arr = np.array(p) / 100.0
        r_arr = np.array(r) / 100.0
        return list((5 * p_arr * r_arr) / (4 * p_arr + r_arr) * 100.0)

    local_ctl_f2 = get_f2(local_ctl_prec, local_ctl_sens)
    local_mos_f2 = get_f2(local_mos_prec, local_mos_sens)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    
    # Colors
    colors = {
        'local_ctl': '#ff7f0e',   # Orange (Control BED)
        'local_mos': '#1f77b4'    # Blue (Mosdepth BED)
    }

    # Subplot 1: Sensitivity
    axs[0, 0].plot(stages, local_ctl_sens, marker='s', color=colors['local_ctl'], linewidth=2.0, linestyle='-.', label='Local Consensus (Control BED)')
    axs[0, 0].plot(stages, local_mos_sens, marker='o', color=colors['local_mos'], linewidth=2.5, linestyle='-', label='Local Consensus (Mosdepth BED)')
    axs[0, 0].set_title('Sensitivity (Recall)', fontsize=12, fontweight='bold')
    axs[0, 0].set_ylabel('Percentage (%)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].set_ylim(50, 105)
    
    # Subplot 2: Precision
    axs[0, 1].plot(stages, local_ctl_prec, marker='s', color=colors['local_ctl'], linewidth=2.0, linestyle='-.', label='Local Consensus (Control BED)')
    axs[0, 1].plot(stages, local_mos_prec, marker='o', color=colors['local_mos'], linewidth=2.5, linestyle='-', label='Local Consensus (Mosdepth BED)')
    axs[0, 1].set_title('Precision (PPV)', fontsize=12, fontweight='bold')
    axs[0, 1].set_ylabel('Percentage (%)')
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)
    axs[0, 1].set_ylim(0, 105)

    # Subplot 3: F2 Score
    axs[1, 0].plot(stages, local_ctl_f2, marker='s', color=colors['local_ctl'], linewidth=2.0, linestyle='-.', label='Local Consensus (Control BED)')
    axs[1, 0].plot(stages, local_mos_f2, marker='o', color=colors['local_mos'], linewidth=2.5, linestyle='-', label='Local Consensus (Mosdepth BED)')
    axs[1, 0].set_title(r'F-beta ($\beta=2$)', fontsize=12, fontweight='bold')
    axs[1, 0].set_ylabel('Percentage (%)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs[1, 0].set_ylim(0, 105)

    # Subplot 4: MCC
    axs[1, 1].plot(stages, local_ctl_mcc, marker='s', color=colors['local_ctl'], linewidth=2.0, linestyle='-.', label='Local Consensus (Control BED)')
    axs[1, 1].plot(stages, local_mos_mcc, marker='o', color=colors['local_mos'], linewidth=2.5, linestyle='-', label='Local Consensus (Mosdepth BED)')
    axs[1, 1].set_title('Matthews Correlation Coefficient (MCC)', fontsize=12, fontweight='bold')
    axs[1, 1].set_ylabel('Percentage (%)')
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)
    axs[1, 1].set_ylim(0, 105)

    for ax in axs.flat:
        ax.set_xticklabels(stages, rotation=20, ha='right', fontsize=9)
        
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=11, frameon=True, shadow=True, bbox_to_anchor=(0.5, 0.01))
    plt.suptitle('Target BED Comparison: Control BED vs. Mosdepth BED (CANOES+XHMM)', fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    
    # Save locations
    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    
    out_path = os.path.join(output_dir, "dragen_vs_local_progression.png")
    plt.savefig(out_path, dpi=300)
    print(f"Saved DRAGEN comparison plot to {out_path}")
    
    import shutil
    shutil.copy(out_path, os.path.join(artifact_dir, "dragen_vs_local_progression.png"))
    plt.close()

if __name__ == '__main__':
    main()
