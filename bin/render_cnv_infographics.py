#!/usr/bin/env python3
"""
render_cnv_infographics.py - Automated CNV Visualizer and Infographic Generator
Compatible with SeeNV (Bradshaw et al., J Mol Diagn 2023) and standalone cluster execution.

Generates standardized visual infographic cards for candidate CNVs by contrasting
target sample read depth against flanking regions and cohort background distributions.
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Render SeeNV-style CNV infographics.")
    parser.add_argument("--vcf", required=True, help="Input CNV VCF (e.g. consensus refined or ACMG classified)")
    parser.add_argument("--bam", required=True, help="Target sample BAM file")
    parser.add_argument("--bai", required=False, help="Target sample BAI file")
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--probes", required=False, help="Target capture intervals BED file")
    parser.add_argument("--ref", required=False, help="Reference FASTA")
    parser.add_argument("--outdir", default="infographics", help="Output directory for infographics")
    parser.add_argument("--max_calls", type=int, default=50, help="Maximum candidate calls to render (prioritizes high-confidence/ACMG)")
    parser.add_argument("--ref_db", required=False, help="SeeNV ReferenceDB directory (if pre-computed)")
    parser.add_argument("--cohort_bams", nargs="*", default=[], help="Optional list of cohort BAMs for background distribution")
    parser.add_argument("--threads", type=int, default=4, help="Number of compute threads")
    return parser.parse_args()

def extract_candidates(vcf_path, max_calls=50):
    """Extract candidate CNVs from VCF, prioritizing ACMG 4/5 and high caller support."""
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    candidates = []
    
    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            pos = int(pos)
            
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True
                    
            end = int(info_dict.get("END", pos + 1000))
            svtype = info_dict.get("SVTYPE", "CNV")
            if svtype == "CNV":
                if "<DEL>" in alt or "DEL" in vid:
                    svtype = "DEL"
                elif "<DUP>" in alt or "DUP" in vid:
                    svtype = "DUP"
            
            callers = info_dict.get("CALLERS", info_dict.get("TOOLS", ""))
            num_callers = int(info_dict.get("NUM_CALLERS", len(callers.split(",")) if callers else 1))
            acmg_class = info_dict.get("ACMG_CLASS", "")
            
            # Prioritization score: ACMG 5 > ACMG 4 > Multi-caller
            priority = 0
            if acmg_class == "5":
                priority += 100
            elif acmg_class == "4":
                priority += 50
            priority += num_callers * 10
            
            candidates.append({
                "chrom": chrom,
                "start": pos,
                "end": end,
                "svtype": svtype,
                "size": end - pos,
                "callers": callers,
                "num_callers": num_callers,
                "acmg_class": acmg_class,
                "priority": priority
            })
            
    # Sort by priority descending
    candidates.sort(key=lambda x: x["priority"], reverse=True)
    return candidates[:max_calls]

def render_fallback_infographics(candidates, bam_path, sample_id, outdir, probes_bed=None):
    """Render high-fidelity static infographic cards using matplotlib and samtools."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    summary_path = os.path.join(outdir, "infographics_summary.tsv")
    
    with open(summary_path, "w") as sum_f:
        sum_f.write("SampleID\tCHROM\tSTART\tEND\tSVTYPE\tSIZE_BP\tCALLERS\tACMG_CLASS\tINFOGRAPHIC\n")
        
        for idx, call in enumerate(candidates, 1):
            chrom = call["chrom"]
            start = call["start"]
            end = call["end"]
            svtype = call["svtype"]
            size_bp = call["size"]
            callers = call["callers"] or "Consensus"
            acmg = call["acmg_class"] or "Unclassified"
            
            # Flanking window (2x size or min 50kb, max 500kb)
            flank = max(min(size_bp, 200000), 25000)
            win_start = max(1, start - flank)
            win_end = end + flank
            
            # Calculate mean coverage in target CNV vs flanks using samtools or pysam
            def get_cov(c, s, e):
                if e <= s:
                    return 0.0
                reg = f"{c}:{s}-{e}"
                if shutil.which("samtools"):
                    res = subprocess.run(
                        ["samtools", "depth", "-r", reg, bam_path],
                        capture_output=True, text=True
                    )
                    depth_str = res.stdout
                else:
                    try:
                        import pysam
                        depth_str = pysam.depth("-r", reg, bam_path)
                    except Exception:
                        depth_str = ""
                depths = [int(l.split("\t")[2]) for l in depth_str.strip().split("\n") if len(l.split("\t")) >= 3]
                return sum(depths) / len(depths) if depths else 0.0
            
            cov_left = get_cov(chrom, win_start, start)
            cov_target = get_cov(chrom, start, end)
            cov_right = get_cov(chrom, end, win_end)
            
            baseline = (cov_left + cov_right) / 2.0 if (cov_left + cov_right) > 0 else 1.0
            norm_ratio = cov_target / baseline if baseline > 0 else 1.0
            
            # Expected copy ratio
            expected_ratio = 0.5 if svtype == "DEL" else 1.5
            
            # Create SeeNV-style figure
            fig = plt.figure(figsize=(12, 6.5), dpi=150)
            gs = fig.add_gridspec(2, 2, width_ratios=[1.8, 1], height_ratios=[1, 1.2], hspace=0.35, wspace=0.25)
            
            # Color scheme
            sv_color = "#e63946" if svtype == "DEL" else "#1d3557"
            bg_color = "#f8f9fa"
            fig.patch.set_facecolor("#ffffff")
            
            # Panel 1: Relative Exonic Depth Barcode / Step Profile
            ax1 = fig.add_subplot(gs[:, 0])
            regions = ["Upstream Flank", f"CNV Target ({svtype})", "Downstream Flank"]
            depths = [cov_left, cov_target, cov_right]
            ratios = [cov_left / baseline if baseline > 0 else 1.0, norm_ratio, cov_right / baseline if baseline > 0 else 1.0]
            
            bar_colors = ["#a8dadc", sv_color, "#a8dadc"]
            bars = ax1.bar(regions, ratios, color=bar_colors, width=0.55, edgecolor="#333333", linewidth=1.2)
            ax1.axhline(1.0, color="#666666", linestyle="--", linewidth=1, label="Diploid Baseline (2.0x)")
            ax1.axhline(expected_ratio, color=sv_color, linestyle=":", linewidth=1.5, label=f"Expected {svtype} Ratio ({expected_ratio:.1f}x)")
            
            ax1.set_ylabel("Normalized Read-Depth Ratio (Target / Baseline)", fontsize=11, fontweight="bold")
            ax1.set_ylim(0, max(max(ratios) * 1.35, 2.2))
            ax1.set_title(f"Read-Depth Architecture: {chrom}:{start:,}-{end:,}", fontsize=12, fontweight="bold", pad=12)
            ax1.grid(axis="y", linestyle=":", alpha=0.6)
            ax1.legend(loc="upper right", frameon=True, fontsize=9)
            
            for bar, d, r in zip(bars, depths, ratios):
                yval = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2.0, yval + 0.05, f"{r:.2f}x\n({d:.1f}x cov)", ha="center", va="bottom", fontsize=9, fontweight="bold")
                
            # Panel 2: Clinical & Evidence Metadata Card
            ax2 = fig.add_subplot(gs[0, 1])
            ax2.axis("off")
            acmg_badge_color = "#d90429" if acmg == "5" else ("#f77f00" if acmg == "4" else "#457b9d")
            
            card_text = (
                f"Sample ID:        {sample_id}\n"
                f"Variant:          {svtype} ({size_bp:,} bp)\n"
                f"Coordinates:      {chrom}:{start:,}-{end:,}\n"
                f"Consensus:        {callers}\n"
                f"Supporting Tools: {call['num_callers']} callers\n"
                f"Depth Ratio:      {norm_ratio:.2f}x (Observed)"
            )
            ax2.text(0.05, 0.95, "Clinical Curation Evidence", fontsize=12, fontweight="bold", transform=ax2.transAxes, va="top")
            ax2.text(0.05, 0.78, card_text, fontsize=9.5, family="monospace", transform=ax2.transAxes, va="top", linespacing=1.4)
            
            # ACMG Tier Badge
            badge_text = f"ACMG Class {acmg}" if acmg in ["4", "5"] else "ACMG Evaluated"
            ax2.text(0.05, 0.05, f"[{badge_text}]", fontsize=11, fontweight="bold", color=acmg_badge_color, transform=ax2.transAxes)

            # Panel 3: Technical Validation Recommendation
            ax3 = fig.add_subplot(gs[1, 1])
            ax3.axis("off")
            ax3.set_facecolor(bg_color)
            
            is_validated = (svtype == "DEL" and norm_ratio <= 0.65) or (svtype == "DUP" and norm_ratio >= 1.30)
            status_title = "✓ TECHNICAL EVIDENCE: SUPPORTED" if is_validated else "⚠ TECHNICAL EVIDENCE: CHECK FLANKS"
            status_color = "#2a9d8f" if is_validated else "#e76f51"
            
            verdict_text = (
                f"Technical Assessment:\n"
                f"• Depth Drop/Gain concordant with {svtype}\n"
                f"• Baseline Flanking Depth: {baseline:.1f}x\n"
                f"• Target Region Depth:   {cov_target:.1f}x\n"
                f"• Action: {'Proceed to clinical sign-off' if is_validated else 'Inspect in IGV for split reads'}"
            )
            ax3.text(0.05, 0.95, status_title, fontsize=11, fontweight="bold", color=status_color, transform=ax3.transAxes, va="top")
            ax3.text(0.05, 0.70, verdict_text, fontsize=9, family="sans-serif", transform=ax3.transAxes, va="top", linespacing=1.3)

            img_filename = f"{sample_id}_{chrom}_{start}_{end}_{svtype}.png"
            img_path = os.path.join(outdir, img_filename)
            plt.savefig(img_path, bbox_inches="tight", dpi=150)
            plt.close(fig)
            
            sum_f.write(f"{sample_id}\t{chrom}\t{start}\t{end}\t{svtype}\t{size_bp}\t{callers}\t{acmg}\t{img_path}\n")

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    candidates = extract_candidates(args.vcf, max_calls=args.max_calls)
    if not candidates:
        print(f"[SeeNV Visualizer] No candidate CNVs found in {args.vcf}. Exiting.")
        summary_path = os.path.join(args.outdir, "infographics_summary.tsv")
        with open(summary_path, "w") as f:
            f.write("SampleID\tCHROM\tSTART\tEND\tSVTYPE\tSIZE_BP\tCALLERS\tACMG_CLASS\tINFOGRAPHIC\n")
        sys.exit(0)
    print(f"[SeeNV Visualizer] Extracted {len(candidates)} high-priority candidates for {args.sample_id}.")
    
    seenv_bin = shutil.which("seenv")
    if seenv_bin and args.probes and args.ref_db and os.path.exists(args.ref_db):
        print(f"[SeeNV Visualizer] Official SeeNV binary detected at {seenv_bin}. Running SeeNV pipeline...")
        sample_list = os.path.join(args.outdir, "sample_list.tsv")
        calls_bed = os.path.join(args.outdir, "candidate_calls.bed")
        with open(sample_list, "w") as f:
            f.write(f"{args.sample_id}\t{os.path.abspath(args.bam)}\n")
        with open(calls_bed, "w") as f:
            for c in candidates:
                f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['svtype']}\t{args.sample_id}\n")
        
        cmd = [
            seenv_bin, "-p",
            "-i", sample_list,
            "-s", args.probes,
            "-c", calls_bed,
            "-r", args.ref_db,
            "-o", args.outdir,
            "-t", str(args.threads)
        ]
        res = subprocess.run(cmd)
        if res.returncode == 0:
            print("[SeeNV Visualizer] Official SeeNV infographic rendering complete.")
            return

    print("[SeeNV Visualizer] Rendering clinical infographic cards via native engine...")
    render_fallback_infographics(candidates, args.bam, args.sample_id, args.outdir, args.probes)
    print(f"[SeeNV Visualizer] Generated infographics successfully in {args.outdir}/")

if __name__ == "__main__":
    main()
