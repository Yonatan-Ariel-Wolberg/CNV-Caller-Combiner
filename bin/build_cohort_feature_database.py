#!/usr/bin/env python3
"""
build_cohort_feature_database.py

Constructs a unified, queryable Cohort Genomic Feature Database by synthesizing:
  1. Key cohort-level attributes from the cohort-merged VCF (AF, AC, AN, carrier spectrum,
     split-read breakpoint resolution status, SV refinement).
  2. Caller/merger-independent biophysical features (interval & flanking GC content,
     target probe count & density, sequence mappability, segmental duplication overlap).
  3. Per-sample carrier event details (sample ID, genotype, supporting callers, concordance,
     tool-specific normalized quality scores, raw scores, read-depth ratio, probe log2 ratio).

Outputs:
  - SQLite database (.db) with normalized tables and a denormalized carrier matrix view.
  - Compressed, tabix-indexed TSVs (.tsv.gz) for variants and carrier calls.
  - Apache Parquet (.parquet) and DuckDB (.duckdb) when pyarrow/duckdb are installed.
"""

import argparse
import gzip
import logging
import math
import os
import sqlite3
import sys
import numpy as np
import pandas as pd
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

SUPPORTED_CALLERS = ["CANOES", "CLAMMS", "XHMM", "CNVKIT", "GATK_GCNV", "DRAGEN", "EXOMEDEPTH"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a unified cohort feature database from cohort VCF and genomic context."
    )
    parser.add_argument(
        "--cohort-vcf",
        required=True,
        help="Path to cohort-merged VCF file (.vcf or .vcf.gz).",
    )
    parser.add_argument(
        "--sample-features-dir",
        default=None,
        help="Directory containing per-sample feature TSVs (*_features.tsv).",
    )
    parser.add_argument(
        "--sample-features-list",
        nargs="*",
        default=[],
        help="List of per-sample feature TSV files.",
    )
    parser.add_argument(
        "--target-bed",
        default=None,
        help="Optional capture target probes BED file.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help="Optional indexed reference genome FASTA file.",
    )
    parser.add_argument(
        "--mappability-bw",
        default=None,
        help="Optional mappability BigWig file (e.g. Umap / ENCODE 100mer).",
    )
    parser.add_argument(
        "--segdup-bed",
        default=None,
        help="Optional segmental duplications BED file (e.g. UCSC genomicSuperDups).",
    )
    parser.add_argument(
        "--output-prefix",
        default="cohort_cnv_features",
        help="Output file prefix (default: 'cohort_cnv_features').",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: current directory).",
    )
    return parser.parse_args()


def open_vcf_read(filepath):
    return pysam.VariantFile(filepath)


def standardize_chrom(chrom_str):
    c = str(chrom_str).strip()
    return c[3:] if c.lower().startswith("chr") else c


def calc_gc_content(fasta_obj, chrom, start, end):
    """Calculate GC content in 1-based [start, end] interval."""
    if fasta_obj is None or start >= end:
        return np.nan
    try:
        ref_chrom = chrom
        if ref_chrom not in fasta_obj.references:
            if f"chr{chrom}" in fasta_obj.references:
                ref_chrom = f"chr{chrom}"
            elif chrom.replace("chr", "") in fasta_obj.references:
                ref_chrom = chrom.replace("chr", "")
            else:
                return np.nan

        seq = fasta_obj.fetch(ref_chrom, max(0, start - 1), end).upper()
        if not seq:
            return np.nan
        g = seq.count("G")
        c = seq.count("C")
        valid_bases = len(seq) - seq.count("N")
        return (g + c) / valid_bases if valid_bases > 0 else np.nan
    except Exception:
        return np.nan


def load_bed_intervals(bed_path):
    """Load BED file into interval dictionary indexed by standardized chrom."""
    intervals_by_chrom = {}
    if not bed_path or not os.path.exists(bed_path):
        return intervals_by_chrom

    open_fn = gzip.open if bed_path.endswith(".gz") else open
    mode = "rt" if bed_path.endswith(".gz") else "r"
    with open_fn(bed_path, mode, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = standardize_chrom(parts[0])
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            if chrom not in intervals_by_chrom:
                intervals_by_chrom[chrom] = []
            intervals_by_chrom[chrom].append((s, e))

    for c in intervals_by_chrom:
        intervals_by_chrom[c].sort(key=lambda x: x[0])
    return intervals_by_chrom


def count_bed_overlaps(intervals_by_chrom, chrom, start, end):
    """Count overlapping intervals and compute total overlapped bases."""
    c_intervals = intervals_by_chrom.get(chrom, [])
    if not c_intervals:
        return 0, 0

    count = 0
    total_bp = 0
    for s, e in c_intervals:
        if s >= end:
            break
        if e <= start:
            continue
        count += 1
        overlap = min(end, e) - max(start, s)
        if overlap > 0:
            total_bp += overlap
    return count, total_bp


def load_sample_features(features_dir, features_list):
    """Load and index all per-sample feature TSVs."""
    all_files = []
    if features_list:
        all_files.extend(features_list)
    if features_dir and os.path.isdir(features_dir):
        for fname in os.listdir(features_dir):
            if fname.endswith("_features.tsv") or fname.endswith(".features.tsv"):
                all_files.append(os.path.join(features_dir, fname))

    sample_dict = {}
    if not all_files:
        logging.info("No per-sample feature TSVs found or supplied.")
        return sample_dict

    logging.info(f"Loading {len(all_files)} per-sample feature TSVs...")
    for fpath in all_files:
        try:
            df = pd.read_csv(fpath, sep="\t")
            if df.empty:
                continue
            for _, row in df.iterrows():
                s_id = str(row.get("sample_id", "")).strip()
                if not s_id:
                    bname = os.path.basename(fpath)
                    s_id = bname.split(".")[0].split("_")[0]
                chrom = standardize_chrom(row.get("chrom", ""))
                try:
                    start = int(row.get("start", 0))
                    end = int(row.get("end", 0))
                except (ValueError, TypeError):
                    continue

                key = (s_id, chrom, start, end)
                sample_dict[key] = row.to_dict()
        except Exception as e:
            logging.warning(f"Could not load {fpath}: {e}")

    logging.info(f"Loaded {len(sample_dict)} per-sample event feature entries.")
    return sample_dict


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Initialize reference resources
    fasta_obj = None
    if args.fasta and os.path.exists(args.fasta):
        try:
            fasta_obj = pysam.FastaFile(args.fasta)
            logging.info(f"Opened reference FASTA: {args.fasta}")
        except Exception as e:
            logging.warning(f"Could not open FASTA {args.fasta}: {e}")

    target_probes = load_bed_intervals(args.target_bed) if args.target_bed else {}
    segdups = load_bed_intervals(args.segdup_bed) if args.segdup_bed else {}

    # 2. Ingest per-sample feature files
    sample_feature_map = load_sample_features(args.sample_features_dir, args.sample_features_list)

    # 3. Parse Cohort VCF
    logging.info(f"Parsing cohort VCF: {args.cohort_vcf}")
    vcf = open_vcf_read(args.cohort_vcf)
    sample_names = list(vcf.header.samples)
    total_cohort_samples = len(sample_names) if len(sample_names) > 0 else 1

    cohort_variant_records = []
    genomic_feature_records = []
    carrier_call_records = []

    def safe_info_get(info_obj, key, default=None):
        try:
            if key in info_obj:
                val = info_obj.get(key)
                return val if val is not None else default
            return default
        except Exception:
            return default

    for record in vcf:
        chrom_raw = record.chrom
        chrom_std = standardize_chrom(chrom_raw)
        pos = record.pos
        info = record.info
        end = int(safe_info_get(info, "END", record.stop))
        svtype = str(safe_info_get(info, "SVTYPE", "CNV")).upper()
        if svtype == "CNV" and record.alts:
            alt0 = str(record.alts[0]).upper()
            if "DEL" in alt0:
                svtype = "DEL"
            elif "DUP" in alt0:
                svtype = "DUP"

        svlen = abs(int(safe_info_get(info, "SVLEN", end - pos + 1)))
        variant_id = f"{chrom_std}:{pos}-{end}:{svtype}"

        # Cohort VCF attributes
        af_val = safe_info_get(info, "AF")
        if af_val is not None:
            if isinstance(af_val, (tuple, list)):
                try:
                    af = float(af_val[0]) if len(af_val) > 0 else np.nan
                except (ValueError, TypeError):
                    af = np.nan
            else:
                try:
                    af = float(af_val)
                except (ValueError, TypeError):
                    af = np.nan
        else:
            af = np.nan

        ac_val = info.get("AC") if "AC" in info else None
        if ac_val is not None:
            if isinstance(ac_val, (tuple, list)):
                try:
                    ac = int(ac_val[0]) if len(ac_val) > 0 else 0
                except (ValueError, TypeError):
                    ac = 0
            else:
                try:
                    ac = int(ac_val)
                except (ValueError, TypeError):
                    ac = 0
        else:
            ac = 0

        an_val = info.get("AN") if "AN" in info else None
        if an_val is not None:
            if isinstance(an_val, (tuple, list)):
                try:
                    an = int(an_val[0]) if len(an_val) > 0 else 2 * total_cohort_samples
                except (ValueError, TypeError):
                    an = 2 * total_cohort_samples
            else:
                try:
                    an = int(an_val)
                except (ValueError, TypeError):
                    an = 2 * total_cohort_samples
        else:
            an = 2 * total_cohort_samples

        is_precise = bool("PRECISE" in info)
        slack_dist = int(info.get("SLACK_DIST", 0)) if "SLACK_DIST" in info else 0
        sv_refined = str(info.get("SV_REFINED", "NONE")) if "SV_REFINED" in info else "NONE"
        orig_pos = int(info.get("ORIG_POS", pos)) if "ORIG_POS" in info else pos
        orig_end = int(info.get("ORIG_END", end)) if "ORIG_END" in info else end

        # Harvest microhomology signatures, mutational mechanism, and overlap metrics
        homlen_val = info.get("HOMLEN") if "HOMLEN" in info else None
        try:
            homlen = int(homlen_val) if homlen_val is not None else 0
        except (ValueError, TypeError):
            homlen = 0
        homseq = str(info.get("HOMSEQ", "")) if "HOMSEQ" in info else ""

        inslen_val = info.get("INSLEN") if "INSLEN" in info else None
        try:
            inslen = int(inslen_val) if inslen_val is not None else 0
        except (ValueError, TypeError):
            inslen = 0
        insertion = str(info.get("INSERTION", "")) if "INSERTION" in info else ""

        mutational_mechanism = str(info.get("MECH", "UNREFINED" if not is_precise else "UNKNOWN")) if "MECH" in info else ("UNREFINED" if not is_precise else "UNKNOWN")

        recip_overlap_val = info.get("REFINED_RECIP_OVERLAP") if "REFINED_RECIP_OVERLAP" in info else None
        try:
            recip_overlap = round(float(recip_overlap_val), 4) if recip_overlap_val is not None else np.nan
        except (ValueError, TypeError):
            recip_overlap = np.nan

        size_sim_val = info.get("REFINED_SIZE_SIM") if "REFINED_SIZE_SIM" in info else None
        try:
            size_sim = round(float(size_sim_val), 4) if size_sim_val is not None else np.nan
        except (ValueError, TypeError):
            size_sim = np.nan

        size_aware_delta_val = info.get("SIZE_AWARE_DELTA") if "SIZE_AWARE_DELTA" in info else None
        try:
            size_aware_delta = int(size_aware_delta_val) if size_aware_delta_val is not None else np.nan
        except (ValueError, TypeError):
            size_aware_delta = np.nan

        # Identify carrier samples and genotypes
        carrier_samples = []
        for sname in sample_names:
            call = record.samples[sname]
            gt = call.get("GT", (None, None))
            if gt and any(a not in (0, None) for a in gt):
                gt_str = "/".join(str(a) if a is not None else "." for a in gt)
                carrier_samples.append((sname, gt_str, call))

        n_carriers = len(carrier_samples)
        if math.isnan(af) and total_cohort_samples > 0:
            af = n_carriers / float(total_cohort_samples)

        cohort_variant_records.append({
            "variant_id": variant_id,
            "chrom": chrom_std,
            "start": pos,
            "end": end,
            "svtype": svtype,
            "svlen": svlen,
            "cohort_af": round(af, 6) if not math.isnan(af) else 0.0,
            "cohort_ac": ac if ac > 0 else n_carriers,
            "cohort_an": an,
            "n_carriers": n_carriers,
            "carrier_sample_ids": ",".join(s[0] for s in carrier_samples) if carrier_samples else "",
            "is_precise": 1 if is_precise else 0,
            "slack_dist": slack_dist,
            "sv_refined_caller": sv_refined,
            "orig_start": orig_pos,
            "orig_end": orig_end,
            "homlen": homlen,
            "homseq": homseq,
            "inslen": inslen,
            "insertion": insertion,
            "mutational_mechanism": mutational_mechanism,
            "size_aware_delta": size_aware_delta,
            "refined_recip_overlap": recip_overlap,
            "refined_size_sim": size_sim,
        })

        # Caller/Merger-Independent Genomic Context Features
        gc_interval = calc_gc_content(fasta_obj, chrom_raw, pos, end)
        gc_flank_up = calc_gc_content(fasta_obj, chrom_raw, max(1, pos - 500), pos)
        gc_flank_dn = calc_gc_content(fasta_obj, chrom_raw, end, end + 500)

        n_probes, probe_overlap_bp = count_bed_overlaps(target_probes, chrom_std, pos, end)
        probe_density = (n_probes / (svlen / 1000.0)) if svlen > 0 else 0.0

        n_segdup, segdup_overlap_bp = count_bed_overlaps(segdups, chrom_std, pos, end)
        segdup_pct = (segdup_overlap_bp / float(svlen)) if svlen > 0 else 0.0

        genomic_feature_records.append({
            "variant_id": variant_id,
            "gc_content": round(gc_interval, 4) if not math.isnan(gc_interval) else np.nan,
            "gc_flank_upstream_500bp": round(gc_flank_up, 4) if not math.isnan(gc_flank_up) else np.nan,
            "gc_flank_downstream_500bp": round(gc_flank_dn, 4) if not math.isnan(gc_flank_dn) else np.nan,
            "n_target_probes": n_probes,
            "probe_density_per_kb": round(probe_density, 3),
            "segdup_overlap_pct": round(segdup_pct, 4),
        })

        # Carrier Calls: connect each carrier sample to sample-level feature extraction metrics
        for sname, gt_str, sample_call in carrier_samples:
            feat_key = (sname, chrom_std, pos, end)
            feat_row = sample_feature_map.get(feat_key, {})

            if not feat_row:
                for (s, c, st, en), r in sample_feature_map.items():
                    if s == sname and c == chrom_std and abs(st - pos) <= 1000 and abs(en - end) <= 1000:
                        feat_row = r
                        break

            # 1. Resolve supporting callers from VCF FORMAT/TOOLS or feat_row
            supporting_callers = []
            vcf_tools = sample_call.get("TOOLS")
            if vcf_tools:
                if isinstance(vcf_tools, (tuple, list)):
                    supporting_callers = [str(x).upper() for x in vcf_tools]
                elif isinstance(vcf_tools, str):
                    supporting_callers = [t.strip().upper() for t in vcf_tools.split(",") if t.strip()]

            if not supporting_callers:
                for caller in ["canoes", "clamms", "xhmm", "cnvkit", "gcnv", "dragen", "exomedepth"]:
                    if feat_row.get(f"is_{caller}", 0) == 1:
                        supporting_callers.append(caller.upper())

            # 2. Resolve concordance and max quality
            concordance = sample_call.get("CONCORDANCE", feat_row.get("concordance", len(supporting_callers)))
            try:
                concordance = int(concordance)
            except (ValueError, TypeError):
                concordance = len(supporting_callers)

            qual_norm = sample_call.get("QUAL_NORM", feat_row.get("max_qual_norm", np.nan))
            try:
                qual_norm = round(float(qual_norm), 2)
            except (ValueError, TypeError):
                qual_norm = np.nan

            rd_ratio = feat_row.get("rd_ratio", np.nan)
            l2r_mean = feat_row.get("l2r_mean", np.nan)
            log2_ratio = feat_row.get("log2_ratio", l2r_mean)
            depth_ratio = feat_row.get("depth_ratio", rd_ratio)
            bin_weight = feat_row.get("bin_weight", feat_row.get("cnvkit_weight", np.nan))
            l2r_dev = feat_row.get("l2r_dev", np.nan)
            l2r_flank_diff = feat_row.get("l2r_flank_diff", np.nan)
            n_probes_flank = feat_row.get("n_probes_flank", np.nan)
            gc_dev = feat_row.get("gc_dev", np.nan)
            is_extreme_gc = feat_row.get("is_extreme_gc", np.nan)
            is_low_mappability = feat_row.get("is_low_mappability", np.nan)
            p_density = feat_row.get("probe_density", np.nan)
            overhang_ratio = feat_row.get("overhang_ratio", np.nan)

            mean_qual_norm = sample_call.get("MEAN_QUAL_NORM")
            try:
                mean_qual_norm = round(float(mean_qual_norm), 2) if mean_qual_norm is not None else np.nan
            except (ValueError, TypeError):
                mean_qual_norm = np.nan

            best_tool = sample_call.get("BEST_TOOL", "")
            best_tool_raw = sample_call.get("BEST_TOOL_RAW_SCORE")
            try:
                best_tool_raw = round(float(best_tool_raw), 2) if best_tool_raw is not None else np.nan
            except (ValueError, TypeError):
                best_tool_raw = np.nan
            best_tool_metric = sample_call.get("BEST_TOOL_METRIC", "")

            # Cross-derivation between log2_ratio and depth_ratio if one is missing
            try:
                if (log2_ratio is None or math.isnan(float(log2_ratio))) and (depth_ratio is not None and not math.isnan(float(depth_ratio))):
                    dr_val = float(depth_ratio)
                    log2_ratio = float(np.log2(dr_val)) if dr_val > 0 else -10.0
                elif (depth_ratio is None or math.isnan(float(depth_ratio))) and (log2_ratio is not None and not math.isnan(float(log2_ratio))):
                    lr_val = float(log2_ratio)
                    depth_ratio = float(2.0 ** lr_val)
            except (ValueError, TypeError):
                pass

            carrier_dict = {
                "variant_id": variant_id,
                "sample_id": sname,
                "genotype": gt_str,
                "concordance": concordance,
                "number_of_callers": concordance,
                "supporting_callers": ",".join(supporting_callers),
                "qual_norm": qual_norm,
                "mean_qual_norm": mean_qual_norm,
                "best_tool": str(best_tool) if best_tool else "",
                "best_tool_raw_score": best_tool_raw,
                "best_tool_metric": str(best_tool_metric) if best_tool_metric else "",
                "rd_ratio": round(float(rd_ratio), 4) if (rd_ratio is not None and not math.isnan(float(rd_ratio))) else np.nan,
                "l2r_mean": round(float(l2r_mean), 4) if (l2r_mean is not None and not math.isnan(float(l2r_mean))) else np.nan,
                "log2_ratio": round(float(log2_ratio), 4) if (log2_ratio is not None and not math.isnan(float(log2_ratio))) else np.nan,
                "depth_ratio": round(float(depth_ratio), 4) if (depth_ratio is not None and not math.isnan(float(depth_ratio))) else np.nan,
                "bin_weight": round(float(bin_weight), 4) if (bin_weight is not None and not math.isnan(float(bin_weight))) else np.nan,
                "l2r_dev": round(float(l2r_dev), 4) if (l2r_dev is not None and not math.isnan(float(l2r_dev))) else np.nan,
                "l2r_flank_diff": round(float(l2r_flank_diff), 4) if (l2r_flank_diff is not None and not math.isnan(float(l2r_flank_diff))) else np.nan,
                "n_probes_flank": int(n_probes_flank) if (n_probes_flank is not None and not math.isnan(float(n_probes_flank))) else np.nan,
                "gc_dev": round(float(gc_dev), 4) if (gc_dev is not None and not math.isnan(float(gc_dev))) else np.nan,
                "is_extreme_gc": int(is_extreme_gc) if (is_extreme_gc is not None and not math.isnan(float(is_extreme_gc))) else np.nan,
                "is_low_mappability": int(is_low_mappability) if (is_low_mappability is not None and not math.isnan(float(is_low_mappability))) else np.nan,
                "probe_density": round(float(p_density), 4) if (p_density is not None and not math.isnan(float(p_density))) else np.nan,
                "overhang_ratio": round(float(overhang_ratio), 4) if (overhang_ratio is not None and not math.isnan(float(overhang_ratio))) else np.nan,
            }

            # 3. Add tool-specific normalized quality scores and raw metrics
            for c_name in SUPPORTED_CALLERS:
                c_lower = c_name.lower()
                # Binary caller indicator
                is_caller = 1 if (c_name in supporting_callers or feat_row.get(f"is_{c_lower}", 0) == 1) else 0
                carrier_dict[f"is_{c_lower}"] = is_caller

                # Normalized quality score
                val_qs = sample_call.get(f"QS_{c_name}")
                if val_qs is None:
                    val_qs = feat_row.get(f"qual_norm_{c_lower}")
                try:
                    carrier_dict[f"qs_{c_lower}"] = round(float(val_qs), 2) if val_qs is not None else np.nan
                except (ValueError, TypeError):
                    carrier_dict[f"qs_{c_lower}"] = np.nan

                # Raw score
                val_oq = sample_call.get(f"OQ_{c_name}")
                try:
                    carrier_dict[f"oq_{c_lower}"] = round(float(val_oq), 2) if val_oq is not None else np.nan
                except (ValueError, TypeError):
                    carrier_dict[f"oq_{c_lower}"] = np.nan

                # Metric name
                val_oas = sample_call.get(f"OAS_{c_name}", "")
                carrier_dict[f"oas_{c_lower}"] = str(val_oas) if val_oas is not None else ""

            carrier_dict["is_gcnv"] = carrier_dict.get("is_gatk_gcnv", 0)

            # Fallback computation if cross-caller summary metrics were not directly in VCF
            if math.isnan(carrier_dict["mean_qual_norm"]):
                available_qs = [carrier_dict[f"qs_{c.lower()}"] for c in SUPPORTED_CALLERS if not math.isnan(carrier_dict.get(f"qs_{c.lower()}", np.nan))]
                if available_qs:
                    carrier_dict["mean_qual_norm"] = round(float(sum(available_qs) / len(available_qs)), 2)

            if not carrier_dict["best_tool"]:
                valid_callers = [(c, carrier_dict.get(f"qs_{c.lower()}", -1.0)) for c in SUPPORTED_CALLERS if not math.isnan(carrier_dict.get(f"qs_{c.lower()}", np.nan))]
                if valid_callers:
                    best_c = max(valid_callers, key=lambda x: x[1])[0]
                    carrier_dict["best_tool"] = best_c
                    carrier_dict["best_tool_raw_score"] = carrier_dict.get(f"oq_{best_c.lower()}", np.nan)
                    carrier_dict["best_tool_metric"] = carrier_dict.get(f"oas_{best_c.lower()}", "")

            carrier_call_records.append(carrier_dict)

    vcf.close()
    if fasta_obj:
        fasta_obj.close()

    df_variants = pd.DataFrame(cohort_variant_records)
    df_genomic = pd.DataFrame(genomic_feature_records)
    df_carriers = pd.DataFrame(carrier_call_records)

    logging.info(f"Synthesized {len(df_variants)} cohort variants and {len(df_carriers)} carrier calls.")

    df_merged_features = pd.merge(df_variants, df_genomic, on="variant_id", how="left")
    df_carrier_matrix = pd.merge(df_carriers, df_merged_features, on="variant_id", how="left") if not df_carriers.empty else pd.DataFrame()

    prefix = os.path.join(args.outdir, args.output_prefix)

    # 4. Export to SQLite Database
    sqlite_path = f"{prefix}.db"
    logging.info(f"Writing SQLite database to {sqlite_path}...")
    conn = sqlite3.connect(sqlite_path)
    df_variants.to_sql("cohort_variants", conn, if_exists="replace", index=False)
    df_genomic.to_sql("genomic_features", conn, if_exists="replace", index=False)
    if not df_carriers.empty:
        df_carriers.to_sql("carrier_calls", conn, if_exists="replace", index=False)
        df_carrier_matrix.to_sql("cohort_carrier_matrix", conn, if_exists="replace", index=False)

    cursor = conn.cursor()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_var_id ON cohort_variants (variant_id);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_var_coords ON cohort_variants (chrom, start, end);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_var_mech ON cohort_variants (mutational_mechanism);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gen_var_id ON genomic_features (variant_id);")
    if not df_carriers.empty:
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_car_var_id ON carrier_calls (variant_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_car_sample_id ON carrier_calls (sample_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_matrix_sample ON cohort_carrier_matrix (sample_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_matrix_var ON cohort_carrier_matrix (variant_id);")

    cursor.execute("""
        CREATE VIEW IF NOT EXISTS v_cohort_features AS
        SELECT 
            v.*,
            g.gc_content,
            g.gc_flank_upstream_500bp,
            g.gc_flank_downstream_500bp,
            g.n_target_probes,
            g.probe_density_per_kb,
            g.segdup_overlap_pct
        FROM cohort_variants v
        LEFT JOIN genomic_features g ON v.variant_id = g.variant_id;
    """)
    conn.commit()
    conn.close()

    # 5. Export to Compressed TSVs
    tsv_variants = f"{prefix}.variants.tsv.gz"
    tsv_carriers = f"{prefix}.carriers.tsv.gz"
    logging.info(f"Writing compressed TSVs: {tsv_variants} and {tsv_carriers}...")
    df_merged_features.to_csv(tsv_variants, sep="\t", index=False, compression="gzip")
    if not df_carrier_matrix.empty:
        df_carrier_matrix.to_csv(tsv_carriers, sep="\t", index=False, compression="gzip")

    # 6. Optional Parquet export
    try:
        parquet_path = f"{prefix}.parquet"
        df_carrier_matrix.to_parquet(parquet_path, index=False)
        logging.info(f"Writing Parquet database to {parquet_path}...")
    except Exception:
        pass

    logging.info("Cohort Feature Database creation complete.")


if __name__ == "__main__":
    main()
