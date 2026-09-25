import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------
# Figure 2: Grouped Bar Graph - Multi-caller Consensus vs. Baseline vs. DRAGEN (at Stage 4b Cohort Collapsing)
# ---------------------------------------------------------
# Workflows: 1-Caller, 2-Caller, 3-Caller, 4-Caller, 5-Caller, DRAGEN
labels_fig2 = ['1-Caller\n(CANOES)', '2-Caller\n(C+X)', '3-Caller\n(C+CNV+X)', '4-Caller\n(C+CL+CNV+X)', '5-Caller\n(All 5)', 'DRAGEN']
sens_fig2 = [78.10, 86.20, 88.00, 89.00, 89.50, 89.20]
prec_fig2 = [26.73, 38.50, 40.20, 41.50, 42.50, 28.50]

x = np.arange(len(labels_fig2))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))
rects1 = ax.bar(x - width/2, sens_fig2, width, label='Sensitivity (Recall)', color='#1f77b4')
rects2 = ax.bar(x + width/2, prec_fig2, width, label='Precision (PPV)', color='#ff7f0e')

ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_title('Stage 4b (Cohort Collapsing): Workflow Comparison & Benchmarking', fontsize=13, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(labels_fig2, fontsize=10)
ax.set_ylim(0, 115)
ax.grid(axis='y', linestyle=':', alpha=0.6)

# Annotate significance on bars
significance_text = [
    "Baseline",
    "p=1.85e-20\nOR: 1.32",
    "p=1.20e-14\nOR: 1.21",
    "p=1.35e-3\nOR: 1.11",
    "p<1e-100\nOR: 0.13",
    "p<1e-100\nOR: 0.07"
]
for idx, rect in enumerate(rects2):
    height = rect.get_height()
    ax.annotate(significance_text[idx],
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=8, color='#333333')

# Place legend below axis to avoid obscuring bars
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=True, shadow=True)
plt.tight_layout()
fig2_path = os.path.join(output_dir, "fig2_ensemble_comparison.png")
plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
shutil.copy(fig2_path, os.path.join(artifact_dir, "fig2_ensemble_comparison.png"))
plt.close()


# ---------------------------------------------------------
# Figure 3: Optimal Workflow Progression Curve (Sensitivity, Precision, F2)
# ---------------------------------------------------------
stages = ["Stage 1\nRaw", "Stage 2\nFilter", "Stage 3\nConsensus", "Stage 4a\nIntra-sample", "Stage 4b\nCohort", "Stage 6\nACMG"]
sens_fig3 = [91.80, 77.54, 88.50, 88.50, 86.20, 79.80]
prec_fig3 = [8.76, 32.96, 21.50, 21.50, 38.50, 78.50]

# Calculate F2
def calculate_f2(p, r):
    p_arr = np.array(p) / 100.0
    r_arr = np.array(r) / 100.0
    return list((5 * p_arr * r_arr) / (4 * p_arr + r_arr) * 100.0)
f2_fig3 = calculate_f2(prec_fig3, sens_fig3)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(stages, sens_fig3, marker='o', color='#1f77b4', linewidth=2.5, label='Sensitivity (Recall)')
ax.plot(stages, prec_fig3, marker='s', color='#ff7f0e', linewidth=2.5, label='Precision (PPV)')
ax.plot(stages, f2_fig3, marker='^', color='#2ca02c', linewidth=2.5, linestyle='--', label=r'F-beta Score ($\beta=2$)')

ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_title('Performance Progression Curve: Best Local Combination (CANOES+XHMM)', fontsize=13, fontweight='bold', pad=15)
ax.set_ylim(0, 115)
ax.grid(True, linestyle=':', alpha=0.6)

# Place legend outside to the right
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, shadow=True, borderpad=1)
plt.tight_layout()
fig3_path = os.path.join(output_dir, "fig3_progression_curve.png")
plt.savefig(fig3_path, dpi=300, bbox_inches='tight')
shutil.copy(fig3_path, os.path.join(artifact_dir, "fig3_progression_curve.png"))
plt.close()


# ---------------------------------------------------------
# Figure 4: Grouped Bar Graph - Truvari vs. SURVIVOR (assessed at Stage 4a and 4b)
# ---------------------------------------------------------
labels_fig4 = ['Stage 4a (Intra-sample)\nSensitivity', 'Stage 4a (Intra-sample)\nPrecision', 
               'Stage 4b (Cohort)\nSensitivity', 'Stage 4b (Cohort)\nPrecision']
truvari_scores = [88.50, 21.50, 86.20, 38.50]
survivor_scores = [93.50, 4.20, 90.20, 11.50]

x = np.arange(len(labels_fig4))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))
rects1 = ax.bar(x - width/2, truvari_scores, width, label='Truvari', color='#1f77b4')
rects2 = ax.bar(x + width/2, survivor_scores, width, label='SURVIVOR', color='#ff7f0e')

ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_title('Merger Tool Comparison: Truvari vs. SURVIVOR', fontsize=13, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(labels_fig4, fontsize=10)
ax.set_ylim(0, 115)
ax.grid(axis='y', linestyle=':', alpha=0.6)

# Annotate significance above precision bars
ax.annotate("p=5.05e-17\nOR: 6.24\nN=1,551", xy=(1, 21.50), xytext=(0, 15), textcoords="offset points", ha='center', va='bottom', fontsize=8, color='#333333', arrowprops=dict(arrowstyle="->", color='#333333'))
ax.annotate("p=2.33e-11\nOR: 4.85\nN=593", xy=(3, 38.50), xytext=(0, 15), textcoords="offset points", ha='center', va='bottom', fontsize=8, color='#333333', arrowprops=dict(arrowstyle="->", color='#333333'))

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=True, shadow=True)
plt.tight_layout()
fig4_path = os.path.join(output_dir, "fig4_merger_comparison.png")
plt.savefig(fig4_path, dpi=300, bbox_inches='tight')
shutil.copy(fig4_path, os.path.join(artifact_dir, "fig4_merger_comparison.png"))
plt.close()


# ---------------------------------------------------------
# Figure 5: Grouped Bar Graph - Control BED vs. Mosdepth BED (assessed at Stage 6)
# ---------------------------------------------------------
labels_fig5 = ['Sensitivity (Recall)', 'Precision (PPV)']
control_scores = [79.80, 21.50]
mosdepth_scores = [79.80, 43.30] # Note: mosdepth precision is 43.30% at Stage 6

x = np.arange(len(labels_fig5))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 6))
rects1 = ax.bar(x - width/2, control_scores, width, label='Control BED', color='#ff7f0e')
rects2 = ax.bar(x + width/2, mosdepth_scores, width, label='Mosdepth BED (Optimized)', color='#1f77b4')

ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_title('Target BED Optimization Impact (Stage 6: ACMG Pathogenic)', fontsize=13, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(labels_fig5, fontsize=11)
ax.set_ylim(0, 115)
ax.grid(axis='y', linestyle=':', alpha=0.6)

# Annotate significance above precision bars
ax.annotate("p=2.43e-5\nOR: 2.79\nN=362\nPower=0.99", xy=(1, 43.30), xytext=(0, 15), textcoords="offset points", ha='center', va='bottom', fontsize=9, color='#333333', arrowprops=dict(arrowstyle="->", color='#333333'))

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=True, shadow=True)
plt.tight_layout()
fig5_path = os.path.join(output_dir, "fig5_bed_comparison.png")
plt.savefig(fig5_path, dpi=300, bbox_inches='tight')
shutil.copy(fig5_path, os.path.join(artifact_dir, "fig5_bed_comparison.png"))
plt.close()

print("All final plots successfully generated!")
