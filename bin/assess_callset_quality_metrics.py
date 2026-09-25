#!/usr/bin/env python3
"""
assess_callset_quality_metrics.py

Extracts and assesses call-level and cohort-level quality metrics across any CNV
callset VCF (single-sample, multi-sample, or cohort-level). Computes:
  1. Variant counts & SVTYPE breakdown (DEL, DUP, other, DEL/DUP ratio)
  2. Size distributions (Min, Q25, Median, Mean, Q75, Max, IQR, total Mb footprint)
  3. Target probe support (Median probes/call, single-probe %, multi-probe % via probes BED)
  4. Quality score distributions (Min, Median, Mean, Max QUAL)
  5. Sample-level distribution (Mean, Median, StdDev calls/sample, outlier samples >200 calls)

Usage:
  python3 assess_callset_quality_metrics.py \
      --vcf-files file1.vcf.gz file2.vcf.gz ... \
      --stage-name "Intra-Sample Inter-Tool Merged" \
      --probes-bed probes.bed \
      --output-summary stage_quality_summary.txt \
      --output-csv stage_quality_metrics.csv
"""

import argparse
import glob
import gzip
import os
import sys
from bisect import bisect_left
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Assess CNV callset quality metrics.")
    parser.add_argument("--vcf-files", nargs="+", required=True, help="Input VCF file(s) or glob pattern.")
    parser.add_argument("--stage-name", required=True, help="Descriptive name of the pipeline callset stage.")
    parser.add_argument("--probes-bed", default=None, help="Target capture probes BED file.")
    parser.add_argument("--output-summary", required=True, help="Output text summary report path.")
    parser.add_argument("--output-csv", default=None, help="Optional detailed per-variant CSV path.")
    parser.add_argument("--output-json", default=None, help="Optional JSON summary path.")
    return parser.parse_args()


def standardize_chrom(c):
    c_str = str(c).strip()
    if c_str.lower().startswith("chr"):
        return c_str[3:]
    return c_str


def load_probes(probes_bed_path):
    if not probes_bed_path or not os.path.exists(probes_bed_path):
        return None
    probes_by_chrom = {}
    try:
        with open(probes_bed_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = standardize_chrom(parts[0])
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                if chrom not in probes_by_chrom:
                    probes_by_chrom[chrom] = {"starts": [], "ends": []}
                probes_by_chrom[chrom]["starts"].append(start)
                probes_by_chrom[chrom]["ends"].append(end)

        for chrom in probes_by_chrom:
            paired = sorted(zip(probes_by_chrom[chrom]["starts"], probes_by_chrom[chrom]["ends"]), key=lambda x: x[0])
            probes_by_chrom[chrom]["starts"] = [p[0] for p in paired]
            probes_by_chrom[chrom]["ends"] = [p[1] for p in paired]
        return probes_by_chrom
    except Exception as e:
        sys.stderr.write(f"Warning: Failed to load probes BED {probes_bed_path}: {e}\n")
        return None


def count_overlapping_probes(probes_by_chrom, chrom, start, end):
    if not probes_by_chrom or chrom not in probes_by_chrom:
        return None
    starts = probes_by_chrom[chrom]["starts"]
    ends = probes_by_chrom[chrom]["ends"]
    n = len(starts)
    if n == 0:
        return 0

    idx = bisect_left(starts, start)
    count = 0
    # Search backwards (up to 50kb buffer for maximum target probe length)
    for i in range(idx - 1, -1, -1):
        if starts[i] >= end:
            continue
        if ends[i] > start:
            count += 1
        if starts[i] < start - 50000:
            break

    # Search forwards
    for i in range(idx, n):
        if starts[i] >= end:
            break
        if ends[i] > start:
            count += 1

    return count


def parse_vcf_records(filepath, probes_by_chrom=None):
    records = []
    opener = gzip.open if filepath.endswith(".gz") else open
    sample_ids = []

    try:
        with opener(filepath, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        parts = line.strip().split("\t")
                        if len(parts) > 9:
                            sample_ids = parts[9:]
                        elif len(parts) == 10:
                            sample_ids = [parts[9]]
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue

                chrom = standardize_chrom(parts[0])
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue

                ref = parts[3]
                alt = parts[4]
                qual_str = parts[5]
                qual = float(qual_str) if qual_str != "." and qual_str != "" else None
                filter_status = parts[6]
                info = parts[7]

                # Determine END coordinate
                end = pos + max(1, len(ref))
                svlen = None
                svtype = "CNV"

                # Parse INFO tags
                for tag in info.split(";"):
                    if tag.startswith("END="):
                        try:
                            end = int(tag.split("=")[1])
                        except ValueError:
                            pass
                    elif tag.startswith("SVLEN="):
                        try:
                            svlen = abs(int(tag.split("=")[1]))
                        except ValueError:
                            pass
                    elif tag.startswith("SVTYPE="):
                        svtype = tag.split("=")[1].strip().upper()

                # Infer SVTYPE from ALT if needed
                if alt in ("<DEL>", "DEL"):
                    svtype = "DEL"
                elif alt in ("<DUP>", "DUP"):
                    svtype = "DUP"
                elif svtype not in ("DEL", "DUP"):
                    if "DEL" in alt:
                        svtype = "DEL"
                    elif "DUP" in alt:
                        svtype = "DUP"

                length = svlen if (svlen is not None and svlen > 0) else max(1, end - pos)

                # Overlapping probes
                probe_count = count_overlapping_probes(probes_by_chrom, chrom, pos, end)

                # Associated samples
                samples_with_variant = []
                if len(parts) > 9 and sample_ids:
                    fmt_keys = parts[8].split(":")
                    gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else -1
                    for s_idx, s_val in enumerate(parts[9:]):
                        if s_idx < len(sample_ids):
                            s_name = sample_ids[s_idx]
                            if gt_idx >= 0:
                                s_gt = s_val.split(":")[gt_idx]
                                if s_gt not in ("0/0", "0|0", "./.", ".", "0"):
                                    samples_with_variant.append(s_name)
                            else:
                                if s_val not in ("./.", ".", ""):
                                    samples_with_variant.append(s_name)
                else:
                    # Single sample derived from filename
                    fname_sample = os.path.basename(filepath).split(".")[0].split("_")[0]
                    samples_with_variant = [fname_sample]

                records.append({
                    "chrom": chrom,
                    "pos": pos,
                    "end": end,
                    "length": length,
                    "svtype": svtype,
                    "qual": qual,
                    "filter": filter_status,
                    "probes": probe_count,
                    "samples": samples_with_variant,
                    "num_samples": len(samples_with_variant)
                })

    except Exception as e:
        sys.stderr.write(f"Warning: Error reading VCF {filepath}: {e}\n")

    return records


def compute_subset_metrics(sub_df, sample_ids_list, total_records):
    count = len(sub_df)
    if count == 0:
        return {
            "count": 0,
            "pct": 0.0,
            "len_min": 0, "len_q25": 0, "len_median": 0, "len_mean": 0.0, "len_q75": 0, "len_max": 0, "len_iqr": 0,
            "footprint_mb": 0.0,
            "probe_median": 0.0, "probe_mean": 0.0, "single_probe_pct": 0.0, "multi_probe_pct": 0.0,
            "qual_min": 0.0, "qual_median": 0.0, "qual_mean": 0.0, "qual_max": 0.0,
            "sample_mean": 0.0, "sample_median": 0.0, "sample_std": 0.0, "outlier_samples": 0
        }

    pct = (count / total_records) * 100.0 if total_records > 0 else 0.0

    lengths = sub_df["length"].values
    len_min = int(np.min(lengths))
    len_q25 = int(np.percentile(lengths, 25))
    len_median = int(np.median(lengths))
    len_mean = float(np.mean(lengths))
    len_q75 = int(np.percentile(lengths, 75))
    len_max = int(np.max(lengths))
    len_iqr = len_q75 - len_q25
    footprint_mb = float(np.sum(lengths)) / 1e6

    probe_vals = [p for p in sub_df["probes"].dropna().values if p is not None]
    if probe_vals:
        probe_median = float(np.median(probe_vals))
        probe_mean = float(np.mean(probe_vals))
        single_probe_pct = (sum(1 for p in probe_vals if p == 1) / len(probe_vals)) * 100.0
        multi_probe_pct = (sum(1 for p in probe_vals if p >= 3) / len(probe_vals)) * 100.0
    else:
        probe_median, probe_mean, single_probe_pct, multi_probe_pct = 0.0, 0.0, 0.0, 0.0

    qual_vals = [q for q in sub_df["qual"].dropna().values if q is not None and not np.isnan(q)]
    if qual_vals:
        qual_min = float(np.min(qual_vals))
        qual_median = float(np.median(qual_vals))
        qual_mean = float(np.mean(qual_vals))
        qual_max = float(np.max(qual_vals))
    else:
        qual_min, qual_median, qual_mean, qual_max = 0.0, 0.0, 0.0, 0.0

    # Sample burden
    sample_sub_counts = {s: 0 for s in sample_ids_list}
    for rec_samples in sub_df["samples"]:
        for s in rec_samples:
            sample_sub_counts[s] = sample_sub_counts.get(s, 0) + 1

    counts_arr = list(sample_sub_counts.values()) if sample_sub_counts else []
    if counts_arr:
        sample_mean = float(np.mean(counts_arr))
        sample_median = float(np.median(counts_arr))
        sample_std = float(np.std(counts_arr))
        outlier_samples = sum(1 for c in counts_arr if c > 200)
    else:
        sample_mean, sample_median, sample_std, outlier_samples = 0.0, 0.0, 0.0, 0

    return {
        "count": count,
        "pct": pct,
        "len_min": len_min, "len_q25": len_q25, "len_median": len_median, "len_mean": len_mean,
        "len_q75": len_q75, "len_max": len_max, "len_iqr": len_iqr,
        "footprint_mb": footprint_mb,
        "probe_median": probe_median, "probe_mean": probe_mean,
        "single_probe_pct": single_probe_pct, "multi_probe_pct": multi_probe_pct,
        "qual_min": qual_min, "qual_median": qual_median, "qual_mean": qual_mean, "qual_max": qual_max,
        "sample_mean": sample_mean, "sample_median": sample_median, "sample_std": sample_std,
        "outlier_samples": outlier_samples
    }


def main():
    args = parse_args()

    # Expand any glob patterns in vcf_files
    all_vcfs = []
    for pattern in args.vcf_files:
        expanded = glob.glob(pattern)
        if expanded:
            all_vcfs.extend(expanded)
        elif os.path.exists(pattern):
            all_vcfs.append(pattern)

    all_vcfs = sorted(list(set(all_vcfs)))
    if not all_vcfs:
        sys.stderr.write(f"Error: No valid VCF files found for pattern: {args.vcf_files}\n")
        with open(args.output_summary, "w", encoding="utf-8") as f:
            f.write(f"# CNV Quality Metrics Assessment: {args.stage_name}\nTotal CNVs: 0\nTotal CNV Calls: 0\nDeletions (DEL): 0 (0.0%)\nDuplications (DUP): 0 (0.0%)\nDEL / DUP Ratio: N/A\nMedian Length: 0 bp\nTotal Genomic Footprint: 0.00 Mb\nMedian Probes per CNV: 0.0\nMean Calls per Sample: 0.0\n")
        if args.output_csv:
            pd.DataFrame(columns=["chrom", "pos", "end", "length", "svtype", "qual", "filter", "probes", "num_samples"]).to_csv(args.output_csv, index=False)
        return

    probes_by_chrom = load_probes(args.probes_bed)

    all_records = []
    sample_call_counts = {}

    for vcf in all_vcfs:
        recs = parse_vcf_records(vcf, probes_by_chrom)
        all_records.extend(recs)
        for r in recs:
            for s in r["samples"]:
                sample_call_counts[s] = sample_call_counts.get(s, 0) + 1

    total_cnvs = len(all_records)
    if total_cnvs == 0:
        with open(args.output_summary, "w", encoding="utf-8") as f:
            f.write(f"# CNV Quality Metrics Assessment: {args.stage_name}\nTotal CNVs: 0\nTotal CNV Calls: 0\nDeletions (DEL): 0 (0.0%)\nDuplications (DUP): 0 (0.0%)\nDEL / DUP Ratio: N/A\nMedian Length: 0 bp\nTotal Genomic Footprint: 0.00 Mb\nMedian Probes per CNV: 0.0\nMean Calls per Sample: 0.0\n")
        if args.output_csv:
            pd.DataFrame(columns=["chrom", "pos", "end", "length", "svtype", "qual", "filter", "probes", "num_samples"]).to_csv(args.output_csv, index=False)
        return

    df = pd.DataFrame(all_records)
    unique_samples = list(sample_call_counts.keys())
    num_unique_samples = len(unique_samples)

    del_df = df[df["svtype"] == "DEL"]
    dup_df = df[df["svtype"] == "DUP"]
    other_df = df[~df["svtype"].isin(["DEL", "DUP"])]

    # Compute subset metrics
    tot_m = compute_subset_metrics(df, unique_samples, total_cnvs)
    del_m = compute_subset_metrics(del_df, unique_samples, total_cnvs)
    dup_m = compute_subset_metrics(dup_df, unique_samples, total_cnvs)
    other_count = len(other_df)
    other_pct = (other_count / total_cnvs) * 100.0 if total_cnvs > 0 else 0.0

    del_dup_ratio = (del_m["count"] / dup_m["count"]) if dup_m["count"] > 0 else float("inf")
    ratio_str = f"{del_dup_ratio:.2f}" if del_dup_ratio != float("inf") else "inf"

    summary_text = f"""====================================================================================================
CNV CALLSET QUALITY METRICS ASSESSMENT REPORT
====================================================================================================
Stage:                          {args.stage_name}
Input Files Evaluated:          {len(all_vcfs)} VCFs
Total CNV Calls:                {total_cnvs:,}
Deletions (DEL):                {del_m['count']:,} ({del_m['pct']:.1f}%)
Duplications (DUP):             {dup_m['count']:,} ({dup_m['pct']:.1f}%)
Other / Complex:                {other_count:,} ({other_pct:.1f}%)
DEL / DUP Ratio:                {ratio_str}
Unique Samples:                 {num_unique_samples:,}

====================================================================================================
3-WAY QUALITY METRIC COMPARISON: TOTAL CNVS vs DELETIONS (DEL) vs DUPLICATIONS (DUP)
====================================================================================================
Metric                          Total CNVs             Deletions (DEL)        Duplications (DUP)    
----------------------------------------------------------------------------------------------------
Total Calls                     {tot_m['count']:<22,d} {f"{del_m['count']:,} ({del_m['pct']:.1f}%)":<22s} {f"{dup_m['count']:,} ({dup_m['pct']:.1f}%)":<22s}
DEL / DUP Ratio                 {ratio_str:<22s} {"-":<22s} {"-":<22s}
Median Size (bp)                {f"{tot_m['len_median']:,} bp":<22s} {f"{del_m['len_median']:,} bp":<22s} {f"{dup_m['len_median']:,} bp":<22s}
Mean Size (bp)                  {f"{tot_m['len_mean']:,.1f} bp":<22s} {f"{del_m['len_mean']:,.1f} bp":<22s} {f"{dup_m['len_mean']:,.1f} bp":<22s}
Interquartile Range (IQR)       {f"{tot_m['len_iqr']:,} bp":<22s} {f"{del_m['len_iqr']:,} bp":<22s} {f"{dup_m['len_iqr']:,} bp":<22s}
Min - Max Size (bp)             {f"{tot_m['len_min']:,} - {tot_m['len_max']:,}":<22s} {f"{del_m['len_min']:,} - {del_m['len_max']:,}":<22s} {f"{dup_m['len_min']:,} - {dup_m['len_max']:,}":<22s}
Total Footprint (Mb)            {f"{tot_m['footprint_mb']:.2f} Mb":<22s} {f"{del_m['footprint_mb']:.2f} Mb":<22s} {f"{dup_m['footprint_mb']:.2f} Mb":<22s}
Median Probes per CNV           {f"{tot_m['probe_median']:.1f}":<22s} {f"{del_m['probe_median']:.1f}":<22s} {f"{dup_m['probe_median']:.1f}":<22s}
Mean Probes per CNV             {f"{tot_m['probe_mean']:.1f}":<22s} {f"{del_m['probe_mean']:.1f}":<22s} {f"{dup_m['probe_mean']:.1f}":<22s}
Single-Probe Calls (%)          {f"{tot_m['single_probe_pct']:.1f}%":<22s} {f"{del_m['single_probe_pct']:.1f}%":<22s} {f"{dup_m['single_probe_pct']:.1f}%":<22s}
Multi-Probe Calls (>=3, %)      {f"{tot_m['multi_probe_pct']:.1f}%":<22s} {f"{del_m['multi_probe_pct']:.1f}%":<22s} {f"{dup_m['multi_probe_pct']:.1f}%":<22s}
Median Quality (QUAL)           {f"{tot_m['qual_median']:.1f}":<22s} {f"{del_m['qual_median']:.1f}":<22s} {f"{dup_m['qual_median']:.1f}":<22s}
Mean Quality (QUAL)             {f"{tot_m['qual_mean']:.1f}":<22s} {f"{del_m['qual_mean']:.1f}":<22s} {f"{dup_m['qual_mean']:.1f}":<22s}
Mean Calls per Sample           {f"{tot_m['sample_mean']:.1f}":<22s} {f"{del_m['sample_mean']:.1f}":<22s} {f"{dup_m['sample_mean']:.1f}":<22s}
Median Calls per Sample         {f"{tot_m['sample_median']:.1f}":<22s} {f"{del_m['sample_median']:.1f}":<22s} {f"{dup_m['sample_median']:.1f}":<22s}
====================================================================================================

--------------------------------------------------------------------------------
1. Variant Classification Breakdown
--------------------------------------------------------------------------------
Deletions (DEL):                {del_m['count']:,} ({del_m['pct']:.1f}%)
Duplications (DUP):             {dup_m['count']:,} ({dup_m['pct']:.1f}%)
Other / Complex:                {other_count:,} ({other_pct:.1f}%)
DEL / DUP Ratio:                {ratio_str}

--------------------------------------------------------------------------------
2. Genomic Size & Footprint Distribution (bp)
--------------------------------------------------------------------------------
Median Length:                  {tot_m['len_median']:,} bp
Median Length (DEL):            {del_m['len_median']:,} bp
Median Length (DUP):            {dup_m['len_median']:,} bp
Mean Length:                    {tot_m['len_mean']:,.1f} bp
Mean Length (DEL):              {del_m['len_mean']:,.1f} bp
Mean Length (DUP):              {dup_m['len_mean']:,.1f} bp
Interquartile Range (IQR):      {tot_m['len_iqr']:,} bp
Total Genomic Footprint:        {tot_m['footprint_mb']:.2f} Mb
Total Genomic Footprint (DEL):  {del_m['footprint_mb']:.2f} Mb
Total Genomic Footprint (DUP):  {dup_m['footprint_mb']:.2f} Mb

--------------------------------------------------------------------------------
3. Target Probe Support Metrics
--------------------------------------------------------------------------------
Median Probes per CNV:          {tot_m['probe_median']:.1f}
Median Probes per CNV (DEL):    {del_m['probe_median']:.1f}
Median Probes per CNV (DUP):    {dup_m['probe_median']:.1f}
Mean Probes per CNV:            {tot_m['probe_mean']:.1f}
Single-Probe Calls (1 probe):   {tot_m['single_probe_pct']:.1f}%
Multi-Probe Calls (>=3 probes): {tot_m['multi_probe_pct']:.1f}%

--------------------------------------------------------------------------------
4. Quality Score Distribution (QUAL)
--------------------------------------------------------------------------------
Minimum QUAL:                   {tot_m['qual_min']:.1f}
Median QUAL:                    {tot_m['qual_median']:.1f}
Mean QUAL:                      {tot_m['qual_mean']:.1f}
Maximum QUAL:                   {tot_m['qual_max']:.1f}

--------------------------------------------------------------------------------
5. Sample Burden & Cohort Uniformity
--------------------------------------------------------------------------------
Mean Calls per Sample:          {tot_m['sample_mean']:.1f}
Mean Calls per Sample (DEL):    {del_m['sample_mean']:.1f}
Mean Calls per Sample (DUP):    {dup_m['sample_mean']:.1f}
Median Calls per Sample:        {tot_m['sample_median']:.1f}
StdDev Calls per Sample:        {tot_m['sample_std']:.1f}
Outlier Samples (>200 calls):   {tot_m['outlier_samples']}
================================================================================
"""

    with open(args.output_summary, "w", encoding="utf-8") as f:
        f.write(summary_text)

    # Save detailed CSV if requested
    if args.output_csv:
        export_df = df[["chrom", "pos", "end", "length", "svtype", "qual", "filter", "probes", "num_samples"]]
        export_df.to_csv(args.output_csv, index=False)

    print(summary_text)


if __name__ == "__main__":
    main()
