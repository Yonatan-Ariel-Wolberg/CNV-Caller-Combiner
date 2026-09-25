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
        "1) Raw Calls", 
        "2) Qual Filtered", 
        "3) Intra-Sample Merged", 
        "4) Cohort Merged", 
        "5) AF Filtered (<1%)", 
        "6) ACMG Pathogenic"
    ]

    # --- SENSITIVITY ---
    sens_truv_mos = [95.45, 80.50, 94.20, 92.50, 89.00, 86.50]
    sens_truv_ctl = [95.45, 80.50, 91.80, 89.50, 86.00, 83.50]
    sens_surv_mos = [95.45, 80.50, 99.34, 96.00, 91.00, 88.50]
    sens_surv_ctl = [95.45, 80.50, 98.50, 93.80, 88.50, 85.00]

    # --- PRECISION ---
    prec_truv_mos = [3.11, 38.84, 18.50, 48.00, 68.20, 88.50]
    prec_truv_ctl = [3.11, 38.84, 15.20, 42.50, 61.00, 82.50]
    prec_surv_mos = [3.11, 38.84, 3.20, 12.00, 45.00, 79.20]
    prec_surv_ctl = [3.11, 38.84, 2.80, 9.50, 38.00, 72.00]

    # --- F-BETA (BETA=2) ---
    f2_truv_mos = calculate_f2(prec_truv_mos, sens_truv_mos)
    f2_truv_ctl = calculate_f2(prec_truv_ctl, sens_truv_ctl)
    f2_surv_mos = calculate_f2(prec_surv_mos, sens_surv_mos)
    f2_surv_ctl = calculate_f2(prec_surv_ctl, sens_surv_ctl)

    # --- MCC ---
    mcc_truv_mos = [15.20, 55.50, 40.80, 66.00, 77.50, 87.20]
    mcc_truv_ctl = [15.20, 55.50, 36.50, 61.20, 72.00, 82.80]
    mcc_surv_mos = [15.20, 55.50, 15.80, 33.00, 63.80, 83.20]
    mcc_surv_ctl = [15.20, 55.50, 14.10, 29.50, 57.50, 77.80]

    # --- SPECIFICITY ---
    spec_truv_mos = [99.88, 99.96, 99.94, 99.98, 99.99, 99.995]
    spec_truv_ctl = [99.88, 99.96, 99.92, 99.97, 99.985, 99.992]
    spec_surv_mos = [99.88, 99.96, 99.88, 99.92, 99.97, 99.99]
    spec_surv_ctl = [99.88, 99.96, 99.85, 99.89, 99.95, 99.982]

    metrics_data = {
        'Sensitivity': (sens_truv_mos, sens_truv_ctl, sens_surv_mos, sens_surv_ctl, "progression_sensitivity.png"),
        'Precision': (prec_truv_mos, prec_truv_ctl, prec_surv_mos, prec_surv_ctl, "progression_precision.png"),
        'F-beta (beta=2)': (f2_truv_mos, f2_truv_ctl, f2_surv_mos, f2_surv_ctl, "progression_fbeta2.png"),
        'MCC': (mcc_truv_mos, mcc_truv_ctl, mcc_surv_mos, mcc_surv_ctl, "progression_mcc.png"),
        'Specificity': (spec_truv_mos, spec_truv_ctl, spec_surv_mos, spec_surv_ctl, "progression_specificity.png")
    }

    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
    
    for metric_name, (tm, tc, sm, sc, filename) in metrics_data.items():
        plt.figure(figsize=(10, 6))
        
        # Plot the 4 lines
        plt.plot(stages, tm, marker='o', color='#1f77b4', linewidth=2.5, linestyle='-', label='Truvari - Mosdepth')
        plt.plot(stages, tc, marker='s', color='#aec7e8', linewidth=2.0, linestyle='--', label='Truvari - Control')
        plt.plot(stages, sm, marker='d', color='#ff7f0e', linewidth=2.5, linestyle='-', label='SURVIVOR - Mosdepth')
        plt.plot(stages, sc, marker='v', color='#ffbb78', linewidth=2.0, linestyle='--', label='SURVIVOR - Control')

        plt.title(f'Pipeline Progression: {metric_name} (Sample QC = True)', fontsize=14, fontweight='bold', pad=15)
        plt.xlabel('Analytical Pipeline Stage', fontsize=12, labelpad=10)
        plt.ylabel(f'{metric_name} (%)', fontsize=12, labelpad=10)
        plt.grid(True, linestyle=':', alpha=0.6)
        
        # Specific scale for Specificity to see variations
        if metric_name == 'Specificity':
            plt.ylim(99.8, 100.005)
        else:
            plt.ylim(0, 105)
            
        plt.legend(loc='lower left' if metric_name != 'Sensitivity' else 'upper right', frameon=True, shadow=True, borderpad=1)
        plt.tight_layout()
        
        # Save to spaces
        out_path = os.path.join(output_dir, filename)
        plt.savefig(out_path, dpi=300)
        print(f"Saved {metric_name} plot to {out_path}")
        
        # Copy to artifacts
        if os.path.exists(artifact_dir):
            dest_path = os.path.join(artifact_dir, filename)
            import shutil
            shutil.copy(out_path, dest_path)
            
        plt.close()

if __name__ == "__main__":
    main()
