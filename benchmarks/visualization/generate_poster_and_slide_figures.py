import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"

metrics = ['Sensitivity', 'Precision', 'F-beta2', 'Specificity', 'MCC']
categories = ['1 Caller', '2 Callers', '3 Callers', '4 Callers', '5 Callers']

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
# Helper to render Concordance Card
# ---------------------------------------------------------
def render_concordance(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    proportions = [91.1, 5.0, 1.9, 1.2, 0.8]
    x_idx = np.arange(len(categories))
    bars = ax.bar(x_idx, proportions, color='#1f77b4', width=0.4, edgecolor='#124f78', zorder=3)
    ax.set_ylabel('Proportion of Total Callset (%)', fontsize=11, fontweight='bold')
    ax.set_title('Caller Consensus Concordance (Stage 2: Filtered)', fontsize=13, fontweight='bold', color='#111827', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0, 110)
    ax.grid(axis='y', linestyle=':', alpha=0.4, zorder=0)

    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    table_data = [
        ["33,526", "1,842", "682", "454", "312"],
        ["91.1%", "5.0%", "1.9%", "1.2%", "0.8%"],
        ["Below Threshold", "Above Threshold", "Above Threshold", "Above Threshold", "Above Threshold"]
    ]
    row_labels = ["Probes Called", "Proportion (%)", "Pipeline Status"]
    table = ax.table(cellText=table_data, rowLabels=row_labels, colLabels=categories, loc='bottom', cellLoc='center')
    style_table(table)

    caption = (
        f"$\\bf{{Figure\\ {fig_num}:}}$ Caller concordance profile (Top) and consensus status table (Bottom) "
        "evaluating probe support tiers across exome calling algorithms."
    )
    fig.text(0.5, -0.15, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# Helper to render Ensemble
# ---------------------------------------------------------
def render_ensemble(fig_num, filename):
    fig, (ax_bar, ax_tbl1, ax_tbl2) = plt.subplots(3, 1, figsize=(9.5, 9.8), gridspec_kw={'height_ratios': [2.2, 1, 1]})
    workflows = ['CANOES', 'XHMM', 'GATK-gCNV', 'CLAMMS', 'CNVkit', 'Ensemble (CANOES+XHMM)']
    sens = [64.41, 57.63, 62.71, 37.29, 66.10, 52.54]
    prec = [38.38, 38.64, 38.54, 44.90, 42.00, 70.45]
    fbeta = [56.38, 52.15, 55.43, 38.60, 58.74, 55.36]
    spec = [99.96, 99.96, 99.96, 99.98, 99.96, 99.98]
    mcc = [0.49, 0.46, 0.48, 0.40, 0.52, 0.61]

    x = np.arange(len(workflows))
    width = 0.15

    rects1 = ax_bar.bar(x - 2*width, sens, width, label='Sensitivity', color='#1f77b4', edgecolor='#124f78', zorder=3)
    rects2 = ax_bar.bar(x - width, prec, width, label='Precision', color='#ff7f0e', edgecolor='#ad5202', zorder=3)
    rects3 = ax_bar.bar(x, fbeta, width, label='F-beta2', color='#2ca02c', edgecolor='#1d6f1d', zorder=3)
    rects4 = ax_bar.bar(x + width, spec, width, label='Specificity', color='#d62728', edgecolor='#8c1617', zorder=3)
    rects5 = ax_bar.bar(x + 2*width, [m*100 for m in mcc], width, label='MCC (x100)', color='#9467bd', edgecolor='#5f3f7d', zorder=3)

    ax_bar.set_ylabel('Performance Metric Value (%)', fontsize=11, fontweight='bold')
    ax_bar.set_title('Performance Benchmark: Single Callers vs. Optimal Ensemble', fontsize=13, fontweight='bold', pad=12)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylim(0, 125)
    ax_bar.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax_bar.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 0.98))

    for rects in [rects1, rects2, rects3, rects4]:
        for rect in rects:
            h = rect.get_height()
            ax_bar.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=7, fontweight='bold')

    # Table 1: Sensitivity
    ax_tbl1.axis('off')
    table_data_1 = [
        ["Reference", "0.68", "0.98", "0.005", "0.91", "0.19"],
        ["Reference", "OR = 1.33", "OR = 1.08", "OR = 3.01", "OR = 0.93", "OR = 1.63"],
        ["Reference", "Power = 0.12", "Power = 0.05", "Power = 0.88", "Power = 0.05", "Power = 0.28"]
    ]
    row_labels_1 = ["Sensitivity (p)", "Odds Ratio (OR)", "Statistical Power"]
    table_1 = ax_tbl1.table(cellText=table_data_1, rowLabels=row_labels_1, colLabels=workflows, loc='center', cellLoc='center')
    table_1.scale(1.0, 1.5)
    table_1.auto_set_font_size(False)
    table_1.set_fontsize(7.5)

    # Table 2: Precision
    ax_tbl2.axis('off')
    table_data_2 = [
        ["Reference", "p < 0.001", "p < 0.001", "p = 0.002", "p < 0.001", "p < 0.001"],
        ["Reference", "OR = 3.82", "OR = 3.78", "OR = 2.91", "OR = 3.25", "Reference"],
        ["Reference", "Power = 0.99", "Power = 0.99", "Power = 0.93", "Power = 0.99", "Reference"]
    ]
    row_labels_2 = ["Precision (p)", "Odds Ratio (OR)", "Statistical Power"]
    table_2 = ax_tbl2.table(cellText=table_data_2, rowLabels=row_labels_2, colLabels=workflows, loc='center', cellLoc='center')
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
        "pairwise statistical significance, odds ratios, and statistical power for sensitivity and precision comparisons."
    )
    fig.text(0.5, 0.03, caption, fontsize=10.5, ha='center', va='center', wrap=True)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.08, hspace=0.15)
    p = os.path.join(output_dir, filename)
    plt.savefig(p, dpi=300, bbox_inches='tight')
    shutil.copy(p, os.path.join(artifact_dir, filename))
    plt.close()

# ---------------------------------------------------------
# Helper to render Merger
# ---------------------------------------------------------
def render_merger(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    x_idx = np.arange(len(metrics))
    width = 0.35
    rects1 = ax.bar(x_idx - width/2, truvari_vals, width, label='Truvari Merging (Stage 3)', color='#1f77b4', edgecolor='#124f78', zorder=3)
    rects2 = ax.bar(x_idx + width/2, survivor_vals, width, label='SURVIVOR Merging (Stage 3)', color='#aec7e8', edgecolor='#6b8ba4', zorder=3)

    ax.set_ylabel('Performance Metric Value (%)', fontsize=11, fontweight='bold')
    ax.set_title('SV Merging Engine Comparison: Truvari vs. SURVIVOR', fontsize=13, fontweight='bold', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0, 130)
    ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 0.98))

    for rect in rects1:
        h = rect.get_height()
        ax.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')
    for rect in rects2:
        h = rect.get_height()
        ax.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')

    table_data = [
        ["52.54%", "70.45%", "55.36%", "99.98%", "0.61"],
        ["52.54%", "55.36%", "53.08%", "99.96%", "0.54"],
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
# Helper to render BED
# ---------------------------------------------------------
def render_bed(fig_num, filename):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    x_idx = np.arange(len(metrics))
    rects1 = ax.bar(x_idx - width/2, control_vals, width, label='Control BED (Stage 1)', color='#2ca02c', edgecolor='#1d6f1d', zorder=3)
    rects2 = ax.bar(x_idx + width/2, mosdepth_vals, width, label='Mosdepth Target BED (Stage 1)', color='#98df8a', edgecolor='#5ca84d', zorder=3)

    ax.set_ylabel('Performance Metric Value (%)', fontsize=11, fontweight='bold')
    ax.set_title('Interval BED Comparison: Control BED vs. Mosdepth target BED', fontsize=13, fontweight='bold', pad=12)
    ax.set_xticks(x_idx)
    ax.set_xticklabels([])
    ax.set_ylim(0, 130)
    ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)
    ax.legend(loc='upper right', framealpha=0.9, bbox_to_anchor=(1.0, 0.98))

    for rect in rects1:
        h = rect.get_height()
        ax.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')
    for rect in rects2:
        h = rect.get_height()
        ax.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, fontweight='bold')

    table_data = [
        ["45.76%", "72.97%", "49.09%", "99.98%", "0.57"],
        ["52.54%", "70.45%", "55.36%", "99.98%", "0.61"],
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
# GENERATE THE TWO SETS OF FIGURES
# ---------------------------------------------------------

# A. Concordance Profile
render_concordance(3, "poster_fig3_concordance.png")
render_concordance(2, "slide_fig2_concordance.png")

# B. Ensemble Comparison
render_ensemble(4, "poster_fig4_ensemble.png")
render_ensemble(3, "slide_fig3_ensemble.png")

# C. Merger Comparison
render_merger(5, "poster_fig5_merger.png")
render_merger(5, "slide_fig5_merger.png")

# D. BED Comparison
render_bed(6, "poster_fig6_bed.png")
render_bed(6, "slide_fig6_bed.png")

# E. Progression Curve (Slide fig 4 only)
fig, ax = plt.subplots(figsize=(8.5, 5.8))
stages = ['Stage 1: Raw Callset', 'Stage 2: Local Consensus', 'Stage 3: Sample QC', 'Stage 4: Control BED', 'Stage 5: Merged (Truvari)']
sensitivity_prog = [64.41, 52.54, 52.54, 45.76, 52.54]
precision_prog = [38.38, 48.44, 55.36, 72.97, 70.45]
ax.plot(stages, sensitivity_prog, marker='o', color='#1f77b4', label='Sensitivity', lw=2.2, markersize=7)
ax.plot(stages, precision_prog, marker='s', color='#ff7f0e', label='Precision', lw=2.2, markersize=7)
ax.set_ylabel('Metric Value (%)', fontsize=11, fontweight='bold')
ax.set_title('Performance Progression Curve Across Pipeline Stages', fontsize=13, fontweight='bold', pad=12)
ax.set_ylim(30, 85)
ax.set_xticklabels(stages, rotation=15, ha='right', fontsize=9)
ax.grid(linestyle=':', alpha=0.6)
ax.legend(loc='upper left', framealpha=0.9)
for i, (s, p) in enumerate(zip(sensitivity_prog, precision_prog)):
    ax.annotate(f'{s:.1f}%', (i, s), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8, fontweight='bold', color='#1f77b4')
    ax.annotate(f'{p:.1f}%', (i, p), xytext=(0, -14), textcoords='offset points', ha='center', fontsize=8, fontweight='bold', color='#ff7f0e')
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

# F. Coverage Tracks (Slide fig 1 only)
fig, (ax1, ax2, ax_gene) = plt.subplots(3, 1, figsize=(8, 6.2), gridspec_kw={'height_ratios': [2, 2, 0.6]}, sharex=True)
np.random.seed(42)
wgs_depth = 30 + np.random.normal(0, 1.2, len(x))
ax1.plot(x, wgs_depth, color='#0288d1', lw=1.5)
ax1.fill_between(x, wgs_depth, color='#e0f2fe', alpha=0.4)
ax1.set_ylabel('Depth (x)', fontsize=9, fontweight='bold')
ax1.set_title('Whole Genome Sequencing (WGS) - Uniform Read Depth', fontsize=10, fontweight='bold', color='#111827', pad=6)
ax1.set_ylim(0, 55)
ax1.grid(axis='y', linestyle=':', alpha=0.4)
es_depth = np.zeros(len(x))
for start, end in exons_coords:
    mask = (x >= start) & (x <= end)
    exon_x = x[mask]
    if start == 15:
        noise = 30 + 18 * np.sin((exon_x - start) * (np.pi / (end - start))) + np.random.normal(0, 2, len(exon_x))
        noise[len(noise)//2] += 12
    elif start == 45:
        noise = 8 + np.random.normal(0, 1.5, len(exon_x))
    else:
        noise = 35 + np.random.normal(0, 2.5, len(exon_x))
    es_depth[mask] = noise
ax2.plot(x, es_depth, color='#d32f2f', lw=1.5)
ax2.fill_between(x, es_depth, color='#ffe2e2', alpha=0.4)
ax2.set_ylabel('Depth (x)', fontsize=9, fontweight='bold')
ax2.set_title('Exome Sequencing (ES) - Capture-Induced Coverage Variance', fontsize=10, fontweight='bold', color='#111827', pad=6)
ax2.set_ylim(0, 60)
ax2.grid(axis='y', linestyle=':', alpha=0.4)
ax2.annotate('Spike (Duplication)', xy=(22.5, 48), xytext=(3, 53), arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8), fontsize=7.5, color='#b91c1c')
ax2.annotate('Drop (Deletion)', xy=(52.5, 8), xytext=(33, 22), arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8), fontsize=7.5, color='#b91c1c')
ax_gene.axis('off')
ax_gene.set_xlim(0, 100)
ax_gene.set_ylim(-1, 1)
ax_gene.plot([0, 100], [0, 0], color='#4b5563', lw=1.5, zorder=1)
for start, end in exons_coords:
    rect = patches.Rectangle((start, -0.3), end - start, 0.6, linewidth=1, edgecolor='#1e3a8a', facecolor='#3b82f6', zorder=2)
    ax_gene.add_patch(rect)
    ax_gene.text((start+end)/2, 0.4, 'Exon', color='#1e3a8a', fontsize=8, fontweight='bold', ha='center')
for start, end in introns_coords:
    ax_gene.text((start+end)/2, -0.7, 'Intron', color='#6b7280', fontsize=8, ha='center')
caption_tracks = (
    r"$\bf{Figure\ 1:}$ Read-depth coverage comparison showing continuous uniform depth in WGS (Top) "
    "vs. targeted, capture-induced noise spikes and drops in ES (Bottom) aligned with gene model exons."
)
fig.text(0.5, 0.02, caption_tracks, fontsize=10.5, ha='center', va='center', wrap=True)
plt.tight_layout()
fig.subplots_adjust(bottom=0.1)
p = os.path.join(output_dir, "slide_fig1_coverage_tracks.png")
plt.savefig(p, dpi=300, bbox_inches='tight')
shutil.copy(p, os.path.join(artifact_dir, "slide_fig1_coverage_tracks.png"))
plt.close()

# G. Duplicate Intro Card as poster_fig1_intro_card.png
shutil.copy(os.path.join(output_dir, "poster_introduction_card.png"), os.path.join(output_dir, "poster_fig1_intro_card.png"))
shutil.copy(os.path.join(output_dir, "poster_introduction_card.png"), os.path.join(artifact_dir, "poster_fig1_intro_card.png"))

# H. Duplicate Master Workflow as poster_fig2_master_workflow.jpg
shutil.copy(os.path.join(output_dir, "phd_entire_workflow_poster.jpg"), os.path.join(output_dir, "poster_fig2_master_workflow.jpg"))
shutil.copy(os.path.join(output_dir, "phd_entire_workflow_poster.jpg"), os.path.join(artifact_dir, "poster_fig2_master_workflow.jpg"))

print("Poster and Slide figures generated successfully!")
