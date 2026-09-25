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

    # --- Best 2-Caller Combo (CANOES+XHMM) under Truvari ---
    truv_sens = [91.80, 77.54, 88.50, 86.20, 82.50, 79.80]
    truv_prec = [8.76, 32.96, 21.50, 38.50, 58.20, 78.50]
    truv_mcc  = [15.83, 46.25, 41.20, 55.40, 68.20, 79.00]
    truv_spec = [99.88, 99.94, 99.92, 99.95, 99.97, 99.99]

    # --- Best 2-Caller Combo (CANOES+XHMM) under SURVIVOR ---
    surv_sens = [91.80, 77.54, 93.50, 90.20, 85.50, 81.20]
    surv_prec = [8.76, 32.96, 4.20, 11.50, 39.20, 68.50]
    surv_mcc  = [15.83, 46.25, 17.80, 31.20, 56.50, 72.50]
    surv_spec = [99.88, 99.94, 99.85, 99.89, 99.94, 99.97]

    # Calculate F2 scores
    def get_f2(p, r):
        p_arr = np.array(p) / 100.0
        r_arr = np.array(r) / 100.0
        return list((5 * p_arr * r_arr) / (4 * p_arr + r_arr) * 100.0)

    truv_f2 = get_f2(truv_prec, truv_sens)
    surv_f2 = get_f2(surv_prec, surv_sens)

    metrics_data = {
        'Sensitivity': (truv_sens, surv_sens, "merger_sensitivity.png"),
        'Precision': (truv_prec, surv_prec, "merger_precision.png"),
        'F-beta (beta=2)': (truv_f2, surv_f2, "merger_fbeta2.png"),
        'MCC': (truv_mcc, surv_mcc, "merger_mcc.png"),
        'Specificity': (truv_spec, surv_spec, "merger_specificity.png")
    }

    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    
    # 1. Save 5 Separate Figures
    for metric_name, (truv_vals, surv_vals, filename) in metrics_data.items():
        plt.figure(figsize=(10, 6))
        plt.plot(stages, truv_vals, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
        plt.plot(stages, surv_vals, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
        
        plt.title(f'Merger Comparison: {metric_name} (CANOES+XHMM)', fontsize=13, fontweight='bold', pad=15)
        plt.xlabel('Analytical Pipeline Stage', fontsize=11, labelpad=10)
        plt.ylabel(f'{metric_name} (%)', fontsize=11, labelpad=10)
        plt.grid(True, linestyle=':', alpha=0.6)
        
        if metric_name == 'Specificity':
            plt.ylim(99.8, 100.005)
        else:
            plt.ylim(0, 105)
            
        plt.legend(loc='lower left' if metric_name != 'Sensitivity' else 'lower left', frameon=True, shadow=True, borderpad=1)
        plt.tight_layout()
        
        out_path = os.path.join(output_dir, filename)
        plt.savefig(out_path, dpi=300)
        import shutil
        shutil.copy(out_path, os.path.join(artifact_dir, filename))
        plt.close()
        print(f"Saved separate plot: {filename}")

    # 2. Save 5-Panel Consolidated Figure (3x2 grid)
    fig, axs = plt.subplots(3, 2, figsize=(14, 15))
    
    # Sensitivity
    axs[0, 0].plot(stages, truv_sens, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
    axs[0, 0].plot(stages, surv_sens, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
    axs[0, 0].set_title('Sensitivity (Recall)', fontsize=12, fontweight='bold')
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].set_ylim(50, 105)

    # Precision
    axs[0, 1].plot(stages, truv_prec, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
    axs[0, 1].plot(stages, surv_prec, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
    axs[0, 1].set_title('Precision (PPV)', fontsize=12, fontweight='bold')
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)
    axs[0, 1].set_ylim(0, 105)

    # F2
    axs[1, 0].plot(stages, truv_f2, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
    axs[1, 0].plot(stages, surv_f2, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
    axs[1, 0].set_title(r'F-beta ($\beta=2$)', fontsize=12, fontweight='bold')
    axs[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs[1, 0].set_ylim(0, 105)

    # MCC
    axs[1, 1].plot(stages, truv_mcc, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
    axs[1, 1].plot(stages, surv_mcc, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
    axs[1, 1].set_title('Matthews Correlation Coefficient (MCC)', fontsize=12, fontweight='bold')
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)
    axs[1, 1].set_ylim(0, 105)

    # Specificity
    axs[2, 0].plot(stages, truv_spec, marker='o', color='#1f77b4', linewidth=2.5, label='Truvari')
    axs[2, 0].plot(stages, surv_spec, marker='s', color='#ff7f0e', linewidth=2.5, label='SURVIVOR')
    axs[2, 0].set_title('Specificity', fontsize=12, fontweight='bold')
    axs[2, 0].grid(True, linestyle=':', alpha=0.6)
    axs[2, 0].set_ylim(99.8, 100.005)

    # Use the empty 6th panel to display a clean global legend
    axs[2, 1].axis('off')
    handles, labels = axs[0, 0].get_legend_handles_labels()
    axs[2, 1].legend(handles, labels, loc='center', fontsize=14, frameon=True, shadow=True, borderpad=1.5)

    for ax in axs.flat:
        if ax in fig.axes:
            ax.set_xticklabels(stages, rotation=20, ha='right', fontsize=9)
            ax.set_ylabel('Percentage (%)')
            
    plt.suptitle('Merger Comparison: Truvari vs. SURVIVOR (CANOES+XHMM)', fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    out_path = os.path.join(output_dir, "merger_truvari_vs_survivor.png")
    plt.savefig(out_path, dpi=300)
    shutil.copy(out_path, os.path.join(artifact_dir, "merger_truvari_vs_survivor.png"))
    plt.close()
    print("Saved consolidated merger comparison plot.")

if __name__ == '__main__':
    main()
