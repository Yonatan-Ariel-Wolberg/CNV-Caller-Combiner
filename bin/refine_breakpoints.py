#!/usr/bin/env python3
"""
refine_breakpoints.py

Post-consensus breakpoint refinement engine for CNV-Caller-Combiner.
Reconciles exon-bounded read-depth consensus calls (from Truvari or SURVIVOR)
with base-pair accurate split-read and assembly structural variants
from DRAGEN-SV (*.sv.vcf.gz) or SvABA (*.svaba.sv.vcf.gz).

Works with pysam when installed, or falls back to pure-Python gzip VCF parsing.
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine CNV consensus breakpoints using orthogonal SV calls."
    )
    parser.add_argument(
        "--consensus-vcf",
        required=True,
        help="Input consensus VCF file (SURVIVOR or Truvari merged).",
    )
    parser.add_argument(
        "--sv-vcf",
        required=True,
        help="Input structural variant VCF file (DRAGEN-SV or SvABA).",
    )
    parser.add_argument(
        "--output-vcf",
        required=True,
        help="Output refined VCF file path.",
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Optional text summary report file path.",
    )
    parser.add_argument(
        "--max-flank",
        type=int,
        default=50000,
        help="Maximum search window into flanking unsequenced introns (bp). Default: 50,000.",
    )
    parser.add_argument(
        "--max-delta",
        type=int,
        default=25000,
        help="Maximum allowed coordinate displacement between depth edge and SV breakpoint (bp). Default: 25,000.",
    )
    parser.add_argument(
        "--focal-threshold",
        type=int,
        default=10000,
        help="Size threshold below which a CNV is treated as focal (bp). Default: 10,000.",
    )
    parser.add_argument(
        "--min-focal-delta",
        type=int,
        default=1000,
        help="Minimum coordinate displacement slack allowed for focal CNVs (bp). Default: 1,000.",
    )
    parser.add_argument(
        "--focal-slack-ratio",
        type=float,
        default=1.0,
        help="Multiplier on focal CNV length to determine allowed coordinate displacement. Default: 1.0.",
    )
    parser.add_argument(
        "--min-recip-overlap",
        type=float,
        default=0.20,
        help="Minimum reciprocal overlap required between consensus CNV and SV for focal events. Default: 0.20.",
    )
    parser.add_argument(
        "--allow-filtered-sv",
        action="store_true",
        help="Allow SV calls with non-PASS filters (useful in targeted exome sequencing where breakends encounter high read depth).",
    )
    parser.add_argument(
        "--engine",
        default=None,
        help="Name of SV engine ('DRAGEN_SV' or 'SVABA'). If omitted, inferred from filename.",
    )
    return parser.parse_args()


def open_vcf_read(filepath):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, "r", encoding="utf-8", errors="replace")


def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info_dict[k] = v
        elif item:
            info_dict[item] = True
    return info_dict


def normalize_svtype(alt, info_dict, id_val=None):
    alt_upper = str(alt).upper()
    if "<DEL" in alt_upper or alt_upper == "DEL":
        return "DEL"
    if "<DUP" in alt_upper or alt_upper == "DUP":
        return "DUP"
    if "<INS" in alt_upper or alt_upper == "INS":
        return "INS"
    if "<INV" in alt_upper or alt_upper == "INV":
        return "INV"

    svtype = info_dict.get("SVTYPE", None)
    if svtype:
        s = str(svtype).upper()
        if "DEL" in s:
            return "DEL"
        if "DUP" in s:
            return "DUP"
        if "INS" in s:
            return "INS"
        if "INV" in s:
            return "INV"
        if s != "CNV":
            return s

    if id_val:
        id_upper = str(id_val).upper()
        if "LOSS" in id_upper or "DEL" in id_upper:
            return "DEL"
        if "GAIN" in id_upper or "DUP" in id_upper:
            return "DUP"

    svlen = info_dict.get("SVLEN")
    if svlen is not None:
        try:
            val = int(svlen)
            if val < 0:
                return "DEL"
            elif val > 0:
                return "DUP"
        except (ValueError, TypeError):
            pass

    return "CNV"


def standardize_chrom(c):
    c_str = str(c).strip()
    return c_str[3:] if c_str.lower().startswith("chr") else c_str


def classify_mutational_mechanism(homlen=0, homseq="", inslen=0, insertion=""):
    """
    Classify structural variant mutational mechanism based on breakpoint junction microhomology
    and non-templated insertion signatures (Carvalho & Lupski 2016; Hastings et al. 2009; Kloosterman et al. 2015).

    Taxonomy:
      - COMPLEX_INS: Non-templated insertion > 10 bp or complex junction rearrangement.
      - NAHR_CANDIDATE: Microhomology >= 25 bp (suggests non-allelic homologous recombination).
      - MMBIR_FOSTES: Microhomology 4-24 bp (suggests microhomology-mediated break-induced replication
        or fork stalling and template switching).
      - NHEJ: 0-3 bp microhomology and insertion <= 10 bp (canonical non-homologous end joining).
      - UNKNOWN: Insufficient breakpoint junction data.
    """
    try:
        homlen = int(homlen or 0)
    except (ValueError, TypeError):
        homlen = 0
    try:
        inslen = int(inslen or 0)
    except (ValueError, TypeError):
        inslen = 0

    if homseq and homlen == 0:
        homlen = len(homseq)
    if insertion and inslen == 0:
        inslen = len(insertion)

    if inslen > 10:
        return "COMPLEX_INS"
    if homlen >= 25:
        return "NAHR_CANDIDATE"
    if 4 <= homlen < 25:
        return "MMBIR_FOSTES"
    if 0 <= homlen <= 3:
        return "NHEJ"
    return "UNKNOWN"


def extract_sv_mechanisms_and_signatures(info):
    """Harvest SvABA and DRAGEN-SV microhomology, insertion, and mechanism fields."""
    homseq_raw = info.get("HOMSEQ", "")
    homseq = "" if (homseq_raw is True or homseq_raw is None or str(homseq_raw) in (".", "None")) else str(homseq_raw)

    homlen_raw = info.get("HOMLEN", 0)
    try:
        homlen = int(homlen_raw)
    except (ValueError, TypeError):
        homlen = len(homseq) if homseq else 0
    if homlen == 0 and homseq:
        homlen = len(homseq)

    ins_raw = info.get("INSERTION", info.get("SVINSSEQ", ""))
    insertion = "" if (ins_raw is True or ins_raw is None or str(ins_raw) in (".", "None")) else str(ins_raw)

    inslen_raw = info.get("INSLEN", 0)
    try:
        inslen = int(inslen_raw)
    except (ValueError, TypeError):
        inslen = len(insertion) if insertion else 0
    if inslen == 0 and insertion:
        inslen = len(insertion)

    mech = classify_mutational_mechanism(homlen=homlen, homseq=homseq, inslen=inslen, insertion=insertion)
    return {
        "homseq": homseq,
        "homlen": homlen,
        "insertion": insertion,
        "inslen": inslen,
        "mech": mech,
    }


def load_sv_records(sv_vcf_path, allow_filtered_sv=False):
    """Load SV records indexed by standardized chromosome for fast interval search."""
    sv_records_by_chrom = {}
    with open_vcf_read(sv_vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            filt = parts[6]
            if not allow_filtered_sv:
                filt_tokens = set(filt.split(";"))
                # In targeted exome sequencing, DRAGEN SV flags breakends with MaxDepth
                # because exome coverage (>50x) exceeds 3x whole-chromosome median (~0.1x).
                # MaxDepth is thus an expected exome coverage flag rather than an SV quality failure.
                allowed_filters = {"PASS", ".", "", "MaxDepth"}
                if not (filt_tokens & allowed_filters or "PASS" in filt):
                    continue

            info = parse_info(parts[7])
            c_std = standardize_chrom(chrom)

            end_val = info.get("END")
            if end_val is not None:
                try:
                    end = int(end_val)
                except (ValueError, TypeError):
                    end = pos + 1
            else:
                m = re.search(r'[\[\]]([^:]+):(\d+)[\[\]]', parts[4])
                if m:
                    mate_c = standardize_chrom(m.group(1))
                    try:
                        mate_p = int(m.group(2))
                        if mate_c == c_std:
                            pos, end = min(pos, mate_p), max(pos, mate_p)
                        else:
                            end = pos + 1
                    except ValueError:
                        end = pos + 1
                else:
                    end = pos + 1

            if end < pos:
                pos, end = end, pos

            svtype = normalize_svtype(parts[4], info, id_val=parts[2])
            sig_dict = extract_sv_mechanisms_and_signatures(info)

            if c_std not in sv_records_by_chrom:
                sv_records_by_chrom[c_std] = []

            sv_records_by_chrom[c_std].append({
                "chrom": chrom,
                "pos": pos,
                "end": end,
                "svtype": svtype,
                "ref": parts[3],
                "alt": parts[4],
                "info": info,
                "homseq": sig_dict["homseq"],
                "homlen": sig_dict["homlen"],
                "insertion": sig_dict["insertion"],
                "inslen": sig_dict["inslen"],
                "mech": sig_dict["mech"],
                "raw_line": line,
            })
    return sv_records_by_chrom


def main():
    args = parse_args()

    if not os.path.exists(args.consensus_vcf):
        sys.exit(f"Error: Consensus VCF not found: {args.consensus_vcf}")
    if not os.path.exists(args.sv_vcf):
        sys.exit(f"Error: SV VCF not found: {args.sv_vcf}")

    engine = args.engine
    if not engine:
        sv_lower = os.path.basename(args.sv_vcf).lower()
        if "dragen" in sv_lower:
            engine = "DRAGEN_SV"
        elif "svaba" in sv_lower:
            engine = "SVABA"
        else:
            engine = "ORTHOGONAL_SV"

    # Pre-index SV calls
    sv_by_chrom = load_sv_records(args.sv_vcf, allow_filtered_sv=args.allow_filtered_sv)

    # Read consensus VCF header and process records
    header_lines = []
    has_precise_def = False
    has_imprecise_def = False
    has_sv_refined_def = False
    has_orig_pos_def = False
    has_orig_end_def = False
    has_slack_dist_def = False
    has_size_aware_delta_def = False
    has_recip_overlap_def = False
    has_size_sim_def = False
    has_homseq_def = False
    has_homlen_def = False
    has_insertion_def = False
    has_inslen_def = False
    has_mech_def = False

    temp_out_uncompressed = args.output_vcf.replace(".gz", "") + ".tmp"
    out_f = open(temp_out_uncompressed, "w", encoding="utf-8")

    total_cnvs = 0
    total_focal = 0
    total_macro = 0
    refined_cnvs = 0
    focal_refined = 0
    macro_refined = 0
    unrefined_cnvs = 0
    del_refined = 0
    dup_refined = 0
    mech_counts = {
        "NHEJ": 0,
        "MMBIR_FOSTES": 0,
        "NAHR_CANDIDATE": 0,
        "COMPLEX_INS": 0,
        "UNKNOWN": 0,
    }

    with open_vcf_read(args.consensus_vcf) as f:
        for line in f:
            if line.startswith("##"):
                if "##INFO=<ID=PRECISE," in line:
                    has_precise_def = True
                elif "##INFO=<ID=IMPRECISE," in line:
                    has_imprecise_def = True
                elif "##INFO=<ID=SV_REFINED," in line:
                    has_sv_refined_def = True
                elif "##INFO=<ID=ORIG_POS," in line:
                    has_orig_pos_def = True
                elif "##INFO=<ID=ORIG_END," in line:
                    has_orig_end_def = True
                elif "##INFO=<ID=SLACK_DIST," in line:
                    has_slack_dist_def = True
                elif "##INFO=<ID=SIZE_AWARE_DELTA," in line:
                    has_size_aware_delta_def = True
                elif "##INFO=<ID=REFINED_RECIP_OVERLAP," in line:
                    has_recip_overlap_def = True
                elif "##INFO=<ID=REFINED_SIZE_SIM," in line:
                    has_size_sim_def = True
                elif "##INFO=<ID=HOMSEQ," in line:
                    has_homseq_def = True
                elif "##INFO=<ID=HOMLEN," in line:
                    has_homlen_def = True
                elif "##INFO=<ID=INSERTION," in line:
                    has_insertion_def = True
                elif "##INFO=<ID=INSLEN," in line:
                    has_inslen_def = True
                elif "##INFO=<ID=MECH," in line:
                    has_mech_def = True
                out_f.write(line)
            elif line.startswith("#CHROM"):
                # Inject precision and refinement definitions before column header
                if not has_precise_def:
                    out_f.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variant has base-pair split-read precision">\n')
                if not has_imprecise_def:
                    out_f.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Coordinates are exon-bound approximations from depth callers">\n')
                if not has_sv_refined_def:
                    out_f.write('##INFO=<ID=SV_REFINED,Number=1,Type=String,Description="Source SV algorithm that refined the breakpoint">\n')
                if not has_orig_pos_def:
                    out_f.write('##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description="Original depth-caller start coordinate (1-based)">\n')
                if not has_orig_end_def:
                    out_f.write('##INFO=<ID=ORIG_END,Number=1,Type=Integer,Description="Original depth-caller end coordinate (1-based)">\n')
                if not has_slack_dist_def:
                    out_f.write('##INFO=<ID=SLACK_DIST,Number=1,Type=Integer,Description="Total boundary delta between depth and SV call (bp)">\n')
                if not has_size_aware_delta_def:
                    out_f.write('##INFO=<ID=SIZE_AWARE_DELTA,Number=1,Type=Integer,Description="Maximum allowed coordinate displacement for variant size (bp)">\n')
                if not has_recip_overlap_def:
                    out_f.write('##INFO=<ID=REFINED_RECIP_OVERLAP,Number=1,Type=Float,Description="Reciprocal overlap between consensus CNV and refining SV">\n')
                if not has_size_sim_def:
                    out_f.write('##INFO=<ID=REFINED_SIZE_SIM,Number=1,Type=Float,Description="Size similarity ratio between consensus CNV and refining SV">\n')
                if not has_homseq_def:
                    out_f.write('##INFO=<ID=HOMSEQ,Number=1,Type=String,Description="Microhomology sequence at structural breakpoint">\n')
                if not has_homlen_def:
                    out_f.write('##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of microhomology at structural breakpoint (bp)">\n')
                if not has_insertion_def:
                    out_f.write('##INFO=<ID=INSERTION,Number=1,Type=String,Description="Non-templated insertion sequence at structural breakpoint">\n')
                if not has_inslen_def:
                    out_f.write('##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Length of non-templated insertion at structural breakpoint (bp)">\n')
                if not has_mech_def:
                    out_f.write('##INFO=<ID=MECH,Number=1,Type=String,Description="Inferred mutational mechanism (NHEJ, MMBIR_FOSTES, NAHR_CANDIDATE, COMPLEX_INS, UNREFINED)">\n')
                out_f.write(line)
            else:
                total_cnvs += 1
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    out_f.write(line)
                    continue

                chrom = parts[0]
                try:
                    c_pos = int(parts[1])
                except ValueError:
                    out_f.write(line)
                    continue

                info = parse_info(parts[7])
                c_end = int(info.get("END", c_pos + 1))
                if c_end < c_pos:
                    c_pos, c_end = c_end, c_pos
                c_len = c_end - c_pos + 1
                c_type = normalize_svtype(parts[4], info, id_val=parts[2])

                is_focal = (c_len < args.focal_threshold)
                if is_focal:
                    total_focal += 1
                    allowed_delta = min(
                        args.max_delta,
                        max(args.min_focal_delta, int(args.focal_slack_ratio * c_len))
                    )
                else:
                    total_macro += 1
                    allowed_delta = args.max_delta

                best_sv = None
                min_delta = float("inf")
                best_ro = 0.0
                best_size_sim = 0.0

                candidates = sv_by_chrom.get(standardize_chrom(chrom), [])
                for sv in candidates:
                    # SVTYPE matching
                    if c_type in ("DEL", "DUP") and sv["svtype"] != c_type:
                        continue

                    # Window matching
                    s_pos = sv["pos"]
                    s_end = sv["end"]
                    s_len = s_end - s_pos + 1

                    if s_pos > c_end + args.max_flank:
                        continue
                    if s_end < c_pos - args.max_flank:
                        continue

                    delta = abs(s_pos - c_pos) + abs(s_end - c_end)
                    if delta > allowed_delta:
                        continue

                    # Overlap and size similarity
                    overlap = max(0, min(c_end, s_end) - max(c_pos, s_pos) + 1)
                    ro = overlap / max(c_len, s_len) if max(c_len, s_len) > 0 else 0.0
                    size_sim = min(c_len, s_len) / max(c_len, s_len) if max(c_len, s_len) > 0 else 0.0

                    if is_focal:
                        # Literature size-aware constraint (GIAB Zook et al. 2020, Truvari English et al. 2022, Carvalho & Lupski 2016):
                        # Focal exonic events must physically overlap the candidate SV and meet either
                        # the reciprocal overlap threshold or complete containment (intronic span).
                        encapsulates = (s_pos <= c_pos and s_end >= c_end)
                        if overlap <= 0:
                            continue
                        if ro < args.min_recip_overlap and not encapsulates:
                            continue

                    if delta < min_delta:
                        min_delta = delta
                        best_sv = sv
                        best_ro = ro
                        best_size_sim = size_sim

                if best_sv is not None and min_delta <= allowed_delta:
                    refined_cnvs += 1
                    if is_focal:
                        focal_refined += 1
                    else:
                        macro_refined += 1

                    if c_type == "DEL":
                        del_refined += 1
                    elif c_type == "DUP":
                        dup_refined += 1

                    mech = best_sv["mech"]
                    mech_counts[mech] = mech_counts.get(mech, 0) + 1

                    # Update POS, REF, ALT, and INFO
                    new_pos = best_sv["pos"]
                    new_end = best_sv["end"]
                    parts[1] = str(new_pos)
                    info["END"] = str(new_end)
                    if "SVLEN" in info:
                        new_len = new_end - new_pos
                        info["SVLEN"] = str(-new_len if c_type == "DEL" else new_len)
                    info["ORIG_POS"] = str(c_pos)
                    info["ORIG_END"] = str(c_end)
                    info["PRECISE"] = True
                    info["SV_REFINED"] = engine
                    info["SLACK_DIST"] = str(int(min_delta))
                    info["SIZE_AWARE_DELTA"] = str(int(allowed_delta))
                    info["REFINED_RECIP_OVERLAP"] = f"{best_ro:.4f}"
                    info["REFINED_SIZE_SIM"] = f"{best_size_sim:.4f}"
                    if best_sv["homlen"] > 0 or best_sv["homseq"]:
                        info["HOMLEN"] = str(best_sv["homlen"])
                    if best_sv["homseq"]:
                        info["HOMSEQ"] = str(best_sv["homseq"])
                    if best_sv["inslen"] > 0 or best_sv["insertion"]:
                        info["INSLEN"] = str(best_sv["inslen"])
                    if best_sv["insertion"]:
                        info["INSERTION"] = str(best_sv["insertion"])
                    info["MECH"] = mech
                    if "IMPRECISE" in info:
                        del info["IMPRECISE"]

                    # Reconstruct INFO field
                    info_items = []
                    for k, v in info.items():
                        if v is True:
                            info_items.append(k)
                        else:
                            info_items.append(f"{k}={v}")
                    parts[7] = ";".join(info_items)
                    out_f.write("\t".join(parts) + "\n")
                else:
                    unrefined_cnvs += 1
                    info["IMPRECISE"] = True
                    if "PRECISE" in info:
                        del info["PRECISE"]
                    info["SIZE_AWARE_DELTA"] = str(int(allowed_delta))
                    info["MECH"] = "UNREFINED"
                    info_items = []
                    for k, v in info.items():
                        if v is True:
                            info_items.append(k)
                        else:
                            info_items.append(f"{k}={v}")
                    parts[7] = ";".join(info_items)
                    out_f.write("\t".join(parts) + "\n")

    out_f.close()

    # Compress and index if needed
    if args.output_vcf.endswith(".gz"):
        sorted_tmp = f"{temp_out_uncompressed}.sorted"
        with open(temp_out_uncompressed, "r", encoding="utf-8") as f_in:
            lines = f_in.readlines()
        headers = [line for line in lines if line.startswith("#")]
        records = [line for line in lines if not line.startswith("#")]
        def sort_key(rec):
            parts = rec.split("\t", 3)
            chrom = parts[0]
            pos = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
            return (chrom, pos)
        records.sort(key=sort_key)
        with open(sorted_tmp, "w", encoding="utf-8") as f_out:
            f_out.writelines(headers)
            f_out.writelines(records)

        has_bgzip = shutil.which("bgzip") is not None
        has_tabix = shutil.which("tabix") is not None
        if has_bgzip and has_tabix:
            with open(args.output_vcf, "wb") as f_out_gz:
                subprocess.run(["bgzip", "-c", sorted_tmp], stdout=f_out_gz, check=True)
            subprocess.run(["tabix", "-p", "vcf", "-f", args.output_vcf], check=True)
            if os.path.exists(sorted_tmp):
                os.remove(sorted_tmp)
            if os.path.exists(temp_out_uncompressed):
                os.remove(temp_out_uncompressed)
        else:
            try:
                import pysam
                pysam.tabix_compress(sorted_tmp, args.output_vcf, force=True)
                pysam.tabix_index(args.output_vcf, preset="vcf", force=True)
                if os.path.exists(sorted_tmp):
                    os.remove(sorted_tmp)
                if os.path.exists(temp_out_uncompressed):
                    os.remove(temp_out_uncompressed)
            except Exception:
                # Pure python gzip fallback
                with open(sorted_tmp if os.path.exists(sorted_tmp) else temp_out_uncompressed, "rb") as f_in:
                    with gzip.open(args.output_vcf, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if os.path.exists(sorted_tmp):
                    os.remove(sorted_tmp)
                if os.path.exists(temp_out_uncompressed):
                    os.remove(temp_out_uncompressed)
    else:
        shutil.move(temp_out_uncompressed, args.output_vcf)

    pct_refined = (refined_cnvs / total_cnvs * 100.0) if total_cnvs > 0 else 0.0
    pct_focal_refined = (focal_refined / total_focal * 100.0) if total_focal > 0 else 0.0
    pct_macro_refined = (macro_refined / total_macro * 100.0) if total_macro > 0 else 0.0

    summary_lines = [
        "==================================================",
        " POST-CONSENSUS BREAKPOINT REFINEMENT SUMMARY",
        "==================================================",
        f"SV Engine:                {engine}",
        f"Total Input CNVs:         {total_cnvs}",
        f"  - Focal (< {args.focal_threshold} bp):  {total_focal}",
        f"  - Macro (>= {args.focal_threshold} bp): {total_macro}",
        f"Refined Base-Pair Exact:  {refined_cnvs} ({pct_refined:.2f}%)",
        f"  - Focal Refined:        {focal_refined} ({pct_focal_refined:.2f}%)",
        f"  - Macro Refined:        {macro_refined} ({pct_macro_refined:.2f}%)",
        f"  - Deletions Refined:    {del_refined}",
        f"  - Duplications Refined: {dup_refined}",
        "Mutational Mechanisms at Junctions:",
        f"  - NHEJ (0-3 bp microhomology):         {mech_counts.get('NHEJ', 0)}",
        f"  - MMBIR/FoSTeS (4-24 bp homology):     {mech_counts.get('MMBIR_FOSTES', 0)}",
        f"  - NAHR candidate (>= 25 bp homology):  {mech_counts.get('NAHR_CANDIDATE', 0)}",
        f"  - Complex Insertion (> 10 bp ins):     {mech_counts.get('COMPLEX_INS', 0)}",
        f"  - Unknown / Untagged SV:               {mech_counts.get('UNKNOWN', 0)}",
        f"Unrefined (Exon-Bounds):  {unrefined_cnvs} ({100.0 - pct_refined:.2f}%)",
        "Preserved Sensitivity:    100.0% (No calls discarded)",
        "==================================================",
    ]
    summary_text = "\n".join(summary_lines)
    print(summary_text)

    if args.summary:
        with open(args.summary, "w") as f:
            f.write(summary_text + "\n")


if __name__ == "__main__":
    main()
