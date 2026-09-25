import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"

metrics = ['Sensitivity', 'Precision', 'F-beta2', 'Specificity', 'MCC']
categories = ['Best 1-Caller\n(CANOES)', 'Best 2-Callers\n(CANOES+XHMM)', 'Best 3-Callers\n(C+X+G)', 'Best 4-Callers\n(C+X+G+CL)', '5-Callers Consensus\n(All Local)', 'DRAGEN Caller\n(No Sample QC)']

truvari_vals = [0.5254, 0.7045, 0.5536, 0.9998, 0.61]
survivor_vals = [0.5254, 0.5536, 0.5308, 0.9996, 0.54]
control_vals = [0.4576, 0.7297, 0.4909, 0.9998, 0.57]
mosdepth_vals = [0.5254, 0.7045, 0.5536, 0.9998, 0.61]
width = 0.35

# Shared table styling helper
def style_table(table, text_size=8.5):
    table.scale(1.0, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(text_size)
    for (row, col), cell in table.get_celld().items():
        cell.set_text_props(fontproperties=matplotlib.font_manager.FontProperties(size=text_size))
        if row == 0 or col == -1:
            cell.set_text_props(fontproperties=matplotlib.font_manager.FontProperties(weight='bold', size=text_size))
            cell.set_facecolor('#f2f2f2')
        else:
            if row in [1, 2]:
                cell.set_text_props(fontproperties=matplotlib.font_manager.FontProperties(size=text_size))
            else:
                cell.set_facecolor('#f8fafc')
        cell.set_edgecolor('#cccccc')

# ---------------------------------------------------------
# A. Concordance Profile (Decimal scale proportions)
# ---------------------------------------------------------
def render_concordance(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    proportions = [0.911, 0.050, 0.019, 0.012, 0.008]
    cats = ['1 Caller', '2 Callers', '3 Callers', '4 Callers', '5 Callers']
    x_idx = np.arange(len(cats))
    bars = ax.bar(x_idx, proportions, color='#1f77b4', width=0.4, edgecolor='#124f78', zorder=3)
    ax.set_ylabel('Proportion of Total Callset', fontsize=11, fontweight='bold')
    ax.set_title('Caller Consensus Concordance (Stage 2: Filtered)', fontsize=13, fontweight='bold', color='#111827', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0.0, 1.1)
    ax.grid(axis='y', linestyle=':', alpha=0.4, zorder=0)

    for bar in bars:
        h = bar.get_height()
        ax.annotate(f'{h:.3f}', xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    table_data = [
        ["33,526", "1,842", "682", "454", "312"],
        ["0.911", "0.050", "0.019", "0.012", "0.008"],
        ["Below Threshold", "Above Threshold", "Above Threshold", "Above Threshold", "Above Threshold"]
    ]
    row_labels = ["Probes Called", "Proportion", "Pipeline Status"]
    table = ax.table(cellText=table_data, rowLabels=row_labels, colLabels=cats, loc='bottom', cellLoc='center')
    style_table(table)

    caption = (
        f"$\\bf{{Figure\\ {fig_num}:}}$ Caller concordance profile (Top) and consensus status table (Bottom) "
        "evaluating probe support tiers across exome calling algorithms (Decimal scale)."
    )
    fig.text(0.5, -0.15, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# B. Ensemble Comparison (Decimal scale)
# ---------------------------------------------------------
def render_tier_ensemble(fig_num, filename):
    fig, (ax_bar, ax_tbl1, ax_tbl2) = plt.subplots(3, 1, figsize=(10.0, 10.2), gridspec_kw={'height_ratios': [2.2, 1, 1]})

    x = np.arange(len(categories))
    w = 0.15

    sens_d = [0.6441, 0.5254, 0.4915, 0.4576, 0.4407, 0.4576]
    prec_d = [0.3838, 0.7045, 0.7250, 0.7297, 0.7222, 0.7105]
    fbeta_d = [0.5638, 0.5536, 0.5234, 0.4909, 0.4738, 0.4902]
    spec_d = [0.9996, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998]
    mcc_d = [0.49, 0.61, 0.59, 0.57, 0.56, 0.57]

    rects1 = ax_bar.bar(x - 2*w, sens_d, w, label='Sensitivity', color='#1f77b4', edgecolor='#124f78', zorder=3)
    rects2 = ax_bar.bar(x - w, prec_d, w, label='Precision', color='#ff7f0e', edgecolor='#ad5202', zorder=3)
    rects3 = ax_bar.bar(x, fbeta_d, w, label='F-beta2', color='#2ca02c', edgecolor='#1d6f1d', zorder=3)
    rects4 = ax_bar.bar(x + w, spec_d, w, label='Specificity', color='#d62728', edgecolor='#8c1617', zorder=3)
    rects5 = ax_bar.bar(x + 2*w, mcc_d, w, label='MCC', color='#9467bd', edgecolor='#5f3f7d', zorder=3)

    ax_bar.set_ylabel('Performance Metric Value', fontsize=11, fontweight='bold')
    ax_bar.set_title('Performance Benchmark: Local Callers Consensus Tiers vs. DRAGEN CNV', fontsize=13, fontweight='bold', pad=12)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylim(0.0, 1.25)
    ax_bar.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax_bar.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 1.15), ncol=5)

    for rects in [rects1, rects2, rects3, rects4, rects5]:
        for rect in rects:
            h = rect.get_height()
            ax_bar.annotate(f'{h:.2f}', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=7, fontweight='bold')

    # Table 1: Sensitivity comparison
    ax_tbl1.axis('off')
    table_data_1 = [
        ["Reference", "p = 0.19", "p = 0.005", "p = 0.005", "p < 0.001", "p = 0.005"],
        ["Reference", "OR = 0.61", "OR = 0.54", "OR = 0.46", "OR = 0.43", "OR = 0.46"],
        ["Reference", "Power = 0.28", "Power = 0.88", "Power = 0.88", "Power = 0.99", "Power = 0.88"]
    ]
    row_labels_1 = ["Sensitivity (p)", "Odds Ratio (OR)", "Statistical Power"]
    table_1 = ax_tbl1.table(cellText=table_data_1, rowLabels=row_labels_1, colLabels=categories, loc='center', cellLoc='center')
    table_1.scale(1.0, 1.5)
    table_1.auto_set_font_size(False)
    table_1.set_fontsize(7.5)

    # Table 2: Precision comparison
    ax_tbl2.axis('off')
    table_data_2 = [
        ["Reference", "p < 0.001", "p < 0.001", "p < 0.001", "p < 0.001", "p < 0.001"],
        ["Reference", "OR = 3.82", "OR = 4.25", "OR = 4.36", "OR = 4.19", "OR = 3.96"],
        ["Reference", "Power = 0.99", "Power = 0.99", "Power = 0.99", "Power = 0.99", "Power = 0.99"]
    ]
    row_labels_2 = ["Precision (p)", "Odds Ratio (OR)", "Statistical Power"]
    table_2 = ax_tbl2.table(cellText=table_data_2, rowLabels=row_labels_2, colLabels=categories, loc='center', cellLoc='center')
    table_2.scale(1.0, 1.5)
    table_2.auto_set_font_size(False)
    table_2.set_fontsize(7.5)

    for t in [table_1, table_2]:
        for (row, col), cell in t.get_celld().items():
            cell.set_text_props(fontproperties=matplotlib.font_manager.FontProperties(size=7.5))
            if row == 0 or col == -1:
                cell.set_text_props(fontproperties=matplotlib.font_manager.FontProperties(weight='bold', size=7.5))
                cell.set_facecolor('#f2f2f2')
            cell.set_edgecolor('#cccccc')

    caption = (
        f"$\\bf{{Figure\\ {fig_num}:}}$ Performance metrics comparison across calling workflows. Lower tables detail "
        "pairwise statistical significance, odds ratios, and statistical power compared to the single-caller baseline (CANOES)."
    )
    fig.text(0.5, 0.02, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.08, hspace=0.18)
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# C. Merger Comparison (Decimal scale)
# ---------------------------------------------------------
def render_merger(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    x_idx = np.arange(len(metrics))
    rects1 = ax.bar(x_idx - width/2, truvari_vals, width, label='Truvari Merging (Stage 3)', color='#1f77b4', edgecolor='#124f78', zorder=3)
    rects2 = ax.bar(x_idx + width/2, survivor_vals, width, label='SURVIVOR Merging (Stage 3)', color='#aec7e8', edgecolor='#6b8ba4', zorder=3)

    ax.set_ylabel('Performance Metric Value', fontsize=11, fontweight='bold')
    ax.set_title('SV Merging Engine Comparison: Truvari vs. SURVIVOR', fontsize=13, fontweight='bold', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0.0, 1.3)
    ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 0.98))

    for rect in rects1:
        h = rect.get_height()
        ax.annotate(f'{h:.3f}', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')
    for rect in rects2:
        h = rect.get_height()
        ax.annotate(f'{h:.3f}', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')

    table_data = [
        ["0.525", "0.705", "0.554", "1.000", "0.61"],
        ["0.525", "0.554", "0.531", "1.000", "0.54"],
        ["p = 1.000", "p = 0.103", "N/A", "p = 0.450", "N/A"],
        ["McNemar (Exact)", "Fisher's Exact", "Descriptor", "Fisher's Exact", "Descriptor"],
        ["OR = 1.00", "OR = 1.92", "N/A", "OR = 2.00", "N/A"]
    ]
    row_labels = ["Truvari Metric", "SURVIVOR Metric", "Significance (p)", "Statistical Test", "Odds Ratio (OR)"]
    table = ax.table(cellText=table_data, rowLabels=row_labels, colLabels=metrics, loc='bottom', cellLoc='center')
    style_table(table)

    caption = (
        f"$\\bf{{Figure\\ {fig_num}:}}$ Comparative evaluation of structural variant merging engines (Truvari vs. SURVIVOR). "
        "Lower table detail McNemar's exact test (Sensitivity) and Fisher's exact test (Precision/Specificity) metrics, odds ratios, and p-values."
    )
    fig.text(0.5, -0.15, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# D. BED Comparison (Decimal scale)
# ---------------------------------------------------------
def render_bed(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    x_idx = np.arange(len(metrics))
    rects1 = ax.bar(x_idx - width/2, control_vals, width, label='Control BED (Stage 1)', color='#2ca02c', edgecolor='#1d6f1d', zorder=3)
    rects2 = ax.bar(x_idx + width/2, mosdepth_vals, width, label='Mosdepth Target BED (Stage 1)', color='#98df8a', edgecolor='#5ca84d', zorder=3)

    ax.set_ylabel('Performance Metric Value', fontsize=11, fontweight='bold')
    ax.set_title('Interval BED Comparison: Control BED vs. Mosdepth target BED', fontsize=13, fontweight='bold', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0.0, 1.3)
    ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 0.98))

    for rect in rects1:
        h = rect.get_height()
        ax.annotate(f'{h:.3f}', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')
    for rect in rects2_b = rects2: # keep variable name robust
        h = rect.get_height()
        ax.annotate(f'{h:.3f}', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')

    table_data = [
        ["0.458", "0.730", "0.491", "1.000", "0.57"],
        ["0.525", "0.705", "0.554", "1.000", "0.61"],
        ["p = 0.125", "p = 1.000", "N/A", "p = 1.000", "N/A"],
        ["McNemar (Exact)", "Fisher's Exact", "Descriptor", "Fisher's Exact", "Descriptor"],
        ["OR = 1.30", "OR = 1.13", "N/A", "OR = 1.00", "N/A"]
    ]
    row_labels = ["Control BED Metric", "Mosdepth Metric", "Significance (p)", "Statistical Test", "Odds Ratio (OR)"]
    table = ax.table(cellText=table_data, rowLabels=row_labels, colLabels=metrics, loc='bottom', cellLoc='center')
    style_table(table)

    caption = (
        f"$\\bf{{Figure\\ {fig_num}:}}$ Performance impact of interval BED generation strategies (Control BED vs. Mosdepth target BED). "
        "Lower table detail McNemar's exact test (Sensitivity) and Fisher's exact test (Precision/Specificity) metrics, odds ratios, and p-values."
    )
    fig.text(0.5, -0.15, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# E. Progression Curve (Slide fig 4 only - Decimal scale)
# ---------------------------------------------------------
def render_progression():
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    stages = ['Stage 1: Raw Callset', 'Stage 2: Local Consensus', 'Stage 3: Sample QC', 'Stage 4: Control BED', 'Stage 5: Merged (Truvari)']
    sensitivity_prog = [0.6441, 0.5254, 0.5254, 0.4576, 0.5254]
    precision_prog = [0.3838, 0.4844, 0.5536, 0.7297, 0.7045]
    ax.plot(stages, sensitivity_prog, marker='o', color='#1f77b4', label='Sensitivity', lw=2.2, markersize=7)
    ax.plot(stages, precision_prog, marker='s', color='#ff7f0e', label='Precision', lw=2.2, markersize=7)
    ax.set_ylabel('Metric Value', fontsize=11, fontweight='bold')
    ax.set_title('Performance Progression Curve Across Pipeline Stages', fontsize=13, fontweight='bold', pad=12)
    ax.set_ylim(0.3, 0.85)
    ax.set_xticklabels(stages, rotation=15, ha='right', fontsize=9)
    ax.grid(linestyle=':', alpha=0.6)
    ax.legend(loc='upper left', framealpha=0.9)
    for i, (s, p) in enumerate(zip(sensitivity_prog, precision_prog)):
        ax.annotate(f'{s:.3f}', (i, s), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8, fontweight='bold', color='#1f77b4')
        ax.annotate(f'{p:.3f}', (i, p), xytext=(0, -14), textcoords='offset points', ha='center', fontsize=8, fontweight='bold', color='#ff7f0e')
    caption_prog = (
        r"$\bf{Figure\ 4:}$ Sequential optimization curves tracking performance progression from single callers "
        "through local consensus filters, sample QC adjustments, control coordinates, and merging engines."
    )
    fig.text(0.5, -0.05, caption_prog, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    p = os.path.join(output_dir, "slide_fig4_progression.png")
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, "slide_fig4_progression.png"))
    plt.close()

# Run all
render_concordance(3, "poster_fig3_concordance.png")
render_concordance(2, "slide_fig2_concordance.png")

render_tier_ensemble(4, "poster_fig4_ensemble.png")
render_tier_ensemble(3, "slide_fig3_ensemble.png")

render_merger(5, "poster_fig5_merger.png")
render_merger(5, "slide_fig5_merger.png")

render_bed(6, "poster_fig6_bed.png")
render_bed(6, "slide_fig6_bed.png")

render_progression()

print("All figures successfully converted to decimal format scale!")
