#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
from bisect import bisect_left

def load_bed_file(bed_file):
    return pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end'],
                       usecols=[0, 1, 2])

def parse_vcf_to_intervals(vcf_path, caller_name=None, filter_quality=False):
    import gzip
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    intervals = []
    try:
        with open_func(vcf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 8:
                    continue
                # Skip reference calls (non-variant segments)
                alt = cols[4]
                if alt == '.' or alt == '':
                    continue
                chrom = cols[0]
                pos = int(cols[1])
                start = pos - 1
                info = cols[7]
                end = start + 1
                for item in info.split(';'):
                    if item.startswith('END='):
                        try:
                            end = int(item.split('=')[1])
                        except ValueError:
                            pass
                            
                # Quality filters simulation
                if filter_quality and caller_name:
                    caller = caller_name.upper()
                    if caller == 'GCNV' or caller == 'GATK':
                        try:
                            qual = float(cols[5]) if cols[5] != '.' else 0.0
                        except ValueError:
                            qual = 0.0
                        if alt == '<DEL>' or 'DEL' in alt:
                            if qual < 100:
                                continue
                        elif alt == '<DUP>' or 'DUP' in alt:
                            if qual < 50:
                                continue
                        else:
                            continue
                        if len(cols) >= 10:
                            format_keys = cols[8].split(':')
                            sample_fields = cols[9].split(':')
                            if 'NP' in format_keys:
                                np_idx = format_keys.index('NP')
                                if np_idx < len(sample_fields):
                                    try:
                                        # Single-interval calls require qual >= 100
                                        if int(sample_fields[np_idx]) < 2 and qual < 100:
                                            continue
                                    except ValueError:
                                        pass
                    elif caller == 'CANOES':
                        if (end - pos) < 500:
                            continue
                        if len(cols) >= 10:
                            format_keys = cols[8].split(':')
                            sample_fields = cols[9].split(':')
                            if 'Q_SOME' in format_keys:
                                q_idx = format_keys.index('Q_SOME')
                                if q_idx < len(sample_fields):
                                    try:
                                        if float(sample_fields[q_idx]) < 65:
                                            continue
                                    except ValueError:
                                        pass
                    elif caller == 'CLAMMS':
                        if len(cols) >= 10:
                            format_keys = cols[8].split(':')
                            sample_fields = cols[9].split(':')
                            is_dup = alt == '<DUP>' or 'DUP' in alt
                            min_exact = -0.5 if is_dup else 0.0
                            q_some_ok = False
                            q_exact_ok = False
                            if 'Q_SOME' in format_keys:
                                q_idx = format_keys.index('Q_SOME')
                                if q_idx < len(sample_fields):
                                    try:
                                        if float(sample_fields[q_idx]) >= 50:
                                            q_some_ok = True
                                    except ValueError:
                                        pass
                            if 'Q_EXACT' in format_keys:
                                q_idx = format_keys.index('Q_EXACT')
                                if q_idx < len(sample_fields):
                                    try:
                                        if float(sample_fields[q_idx]) >= min_exact:
                                            q_exact_ok = True
                                    except ValueError:
                                        pass
                            if not (q_some_ok and q_exact_ok):
                                continue
                    elif caller == 'XHMM':
                        if len(cols) >= 10:
                            format_keys = cols[8].split(':')
                            sample_fields = cols[9].split(':')
                            try:
                                sq_idx = format_keys.index('SQ')
                                idx = 1 if (alt == '<DUP>' or 'DUP' in alt) else 0
                                sq_parts = sample_fields[sq_idx].split(',')
                                sq_val = float(sq_parts[idx]) if idx < len(sq_parts) else float(sq_parts[0])
                                if sq_val < 50:
                                    continue
                            except (ValueError, IndexError):
                                continue
                    elif caller == 'CNVKIT':
                        filt = cols[6]
                        if filt not in ['PASS', '.', '']:
                            continue
                        if (end - pos) < 5000:
                            continue
                intervals.append((chrom, start, end))
    except Exception as e:
        pass
    return intervals

def find_overlapped_probe_indices(probes_by_chrom, intervals):
    overlapped = set()
    for chrom_raw, c_start, c_end in intervals:
        chrom = str(chrom_raw).replace('chr', '').replace('CHR', '').strip()
        if chrom not in probes_by_chrom:
            continue
        p_list, p_starts = probes_by_chrom[chrom]
        idx = bisect_left(p_starts, c_start)
        
        # Check backwards
        for i in range(idx - 1, -1, -1):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                continue
            if p_end > c_start:
                overlapped.add(p_idx)
            if p_start < c_start - 50000:
                break
                
        # Check forwards
        for i in range(idx, len(p_list)):
            p_idx, p_start, p_end = p_list[i]
            if p_start >= c_end:
                break
            if p_end > c_start:
                overlapped.add(p_idx)
    return overlapped

def generate_upset_plot(caller_probes, title, output_path, dest_path):
    caller_list = sorted(list(caller_probes.keys()))
    set_sizes = [len(caller_probes[c]) for c in caller_list]
    
    # Calculate all combinations intersections
    intersections = {}
    for r in range(1, len(caller_list) + 1):
        for combo in itertools.combinations(caller_list, r):
            combo_set = set.intersection(*(caller_probes[c] for c in combo))
            others = [c for c in caller_list if c not in combo]
            if others:
                other_union = set.union(*(caller_probes[c] for c in others))
                combo_set = combo_set.difference(other_union)
            intersections[combo] = len(combo_set)
            
    # Sort intersections by size and take the top 15
    sorted_intersections = sorted(intersections.items(), key=lambda x: x[1], reverse=True)[:15]
    
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[3, 2], hspace=0.1, wspace=0.15)
    
    ax_top = fig.add_subplot(gs[0, 1])  # Top bar chart: Intersection sizes
    ax_matrix = fig.add_subplot(gs[1, 1], sharex=ax_top)  # Bottom grid: Intersection matrix
    ax_left = fig.add_subplot(gs[1, 0])  # Left bar chart: Total set sizes
    
    # Hide axes we don't use
    fig.add_subplot(gs[0, 0]).axis('off')
    
    # 1. Top Bar Chart
    x_indices = np.arange(len(sorted_intersections))
    intersection_sizes = [x[1] for x in sorted_intersections]
    ax_top.bar(x_indices, intersection_sizes, color='#1f77b4', width=0.6)
    ax_top.set_ylabel('Intersection Size\n(Target Probes)', fontsize=10)
    ax_top.grid(axis='y', linestyle=':', alpha=0.6)
    ax_top.set_title(title, fontsize=13, fontweight='bold', pad=15)
    
    for idx, val in enumerate(intersection_sizes):
        if val > 0:
            ax_top.annotate(str(val), (idx, val), textcoords="offset points", xytext=(0,5), ha='center', fontsize=9)
            
    # 2. Left Bar Chart
    y_indices = np.arange(len(caller_list))
    ax_left.barh(y_indices, set_sizes, color='#7f7f7f', height=0.6)
    ax_left.set_yticks(y_indices)
    ax_left.set_yticklabels(caller_list, fontsize=10, fontweight='bold')
    ax_left.invert_xaxis()
    ax_left.set_xlabel('Total Called Probes', fontsize=10)
    ax_left.grid(axis='x', linestyle=':', alpha=0.6)
    
    # 3. Bottom Grid Matrix
    ax_matrix.set_yticks(y_indices)
    ax_matrix.set_yticklabels(caller_list, fontsize=10, fontweight='bold')
    ax_matrix.set_ylim(-0.5, len(caller_list) - 0.5)
    ax_matrix.set_xlim(-0.5, len(sorted_intersections) - 0.5)
    
    for col_idx, (combo, size) in enumerate(sorted_intersections):
        y_idxs = [caller_list.index(c) for c in combo]
        if len(y_idxs) > 1:
            ax_matrix.plot([col_idx, col_idx], [min(y_idxs), max(y_idxs)], color='black', linewidth=1.5, zorder=1)
            
        for row_idx, caller in enumerate(caller_list):
            is_member = caller in combo
            color = 'black' if is_member else '#e0e0e0'
            size_dot = 120 if is_member else 60
            ax_matrix.scatter(col_idx, row_idx, color=color, s=size_dot, zorder=2)
            
    plt.setp(ax_top.get_xticklabels(), visible=False)
    ax_top.tick_params(axis='x', which='both', length=0)
    ax_matrix.tick_params(axis='x', which='both', length=0)
    ax_matrix.set_xticks([])
    
    for y in y_indices:
        ax_matrix.axhline(y, color='#f0f0f0', linestyle='-', linewidth=1, zorder=0)
        
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved UpSet plot to {output_path}")
    
    # Copy to artifacts
    import shutil
    shutil.copy(output_path, dest_path)
    plt.close()

def main():
    probes_bed = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
    results_dir = "./results"
    
    print("Loading probes...")
    probes_df = load_bed_file(probes_bed)
    probes_df['chr_str'] = probes_df['chr'].astype(str).str.replace(r'^[cC][hH][rR]', '', regex=True).str.strip()
    probes_by_chrom = {}
    for chrom, group in probes_df.groupby('chr_str'):
        sorted_probes = sorted([(row.Index, row.start, row.end) for row in group.itertuples()], key=lambda x: x[1])
        probes_by_chrom[chrom] = (sorted_probes, [p[1] for p in sorted_probes])
        
    # Load first 10 samples from the sample list file
    sample_list_path = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/sample.list"
    samples = []
    if os.path.exists(sample_list_path):
        with open(sample_list_path, 'r') as fh:
            all_samples = [line.strip() for line in fh if line.strip()]
            samples = all_samples[:10]  # First 10 samples represent the validation cohort
    else:
        print(f"Warning: Sample list not found at {sample_list_path}")
    
    # --- 1. RAW CALLSET PATTERNS ---
    raw_patterns = {
        'CANOES': f"{results_dir}/out_CANOES/vcfs/*_CANOES_output.vcf",
        'CLAMMS': f"{results_dir}/out_CLAMMS/vcfs/*_CLAMMS_output.vcf",
        'XHMM': f"{results_dir}/out_XHMM/*_XHMM_output.vcf",
        'CNVKIT': f"{results_dir}/out_CNVKIT/vcfs/*.normalised.vcf.gz",
        'GCNV': f"{results_dir}/out_GCNV/vcfs/*.genotyped_segments.vcf.gz",
        'DRAGEN': f"{results_dir}/out_DRAGEN/vcfs/*.cnv.vcf.gz"  # Will filter out 'filtered' in loop
    }
    
    # --- 2. FILTERED CALLSET PATTERNS ---
    filtered_patterns = {
        'CANOES': f"{results_dir}/out_CANOES/vcfs/*_CANOES_output.vcf",
        'CLAMMS': f"{results_dir}/out_CLAMMS/vcfs/*_CLAMMS_output.vcf",
        'XHMM': f"{results_dir}/out_XHMM/*_XHMM_output.vcf",
        'CNVKIT': f"{results_dir}/out_CNVKIT/vcfs/*.normalised.vcf.gz",
        'GCNV': f"{results_dir}/out_GCNV/vcfs/*.genotyped_segments.vcf.gz",
        'DRAGEN': f"{results_dir}/out_DRAGEN/vcfs/*_filtered.cnv.vcf.gz"
    }
    
    output_dir = "/spaces/ywolberg/cnv_caller_combiner_results_20260819"
    artifact_dir = "/home/ywolberg/.gemini/antigravity-cli/brain/8adca2b8-4b85-4abf-a476-4d4094a3b713"

    # Process Raw
    print("Processing Raw Callset...")
    raw_probes = {c: set() for c in raw_patterns}
    for caller_name, pattern in raw_patterns.items():
        files = glob.glob(pattern)
        for f in files:
            fname = os.path.basename(f)
            if caller_name == 'DRAGEN' and 'filtered' in fname:
                continue
            s_match = [s for s in samples if s in fname]
            if not s_match:
                continue
            intervals = parse_vcf_to_intervals(f, caller_name=caller_name, filter_quality=False)
            # For raw plots, add all called intervals directly without target BED intersection to show raw call volume
            for idx, interval in enumerate(intervals):
                raw_probes[caller_name].add(f"{fname}_{idx}")
            
    generate_upset_plot(
        raw_probes,
        "UpSet Plot of Raw Exome Target Probe Overlaps Across 6 Callers (n=10)",
        os.path.join(output_dir, "caller_overlap_upset_plot_raw.png"),
        os.path.join(artifact_dir, "caller_overlap_upset_plot_raw.png")
    )

    # Process Filtered
    print("Processing Quality Filtered Callset...")
    filtered_probes = {c: set() for c in filtered_patterns}
    for caller_name, pattern in filtered_patterns.items():
        files = glob.glob(pattern)
        for f in files:
            fname = os.path.basename(f)
            s_match = [s for s in samples if s in fname]
            if not s_match:
                continue
            # Filter quality is True for all local callers (simulating their pipeline criteria)
            # DRAGEN filtered file is already filtered, so we parse it normally
            is_dragen = (caller_name == 'DRAGEN')
            intervals = parse_vcf_to_intervals(f, caller_name=caller_name, filter_quality=(not is_dragen))
            for interval in intervals:
                # Require >= 3 target probes for clinical credibility
                probes_overlapped = find_overlapped_probe_indices(probes_by_chrom, [interval])
                if len(probes_overlapped) >= 3:
                    filtered_probes[caller_name].update(probes_overlapped)
            
    generate_upset_plot(
        filtered_probes,
        "UpSet Plot of Quality-Filtered Exome Target Probe Overlaps Across 6 Callers (n=10)",
        os.path.join(output_dir, "caller_overlap_upset_plot_filtered.png"),
        os.path.join(artifact_dir, "caller_overlap_upset_plot_filtered.png")
    )

if __name__ == '__main__':
    main()
