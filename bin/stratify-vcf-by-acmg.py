#!/usr/bin/env python3
"""
stratify-vcf-by-acmg.py

Stratifies an AnnotSV-annotated cohort VCF into individual ACMG pathogenicity-class
VCFs (ACMG1, ACMG2, ACMG3, ACMG4, ACMG5, ACMGNA, and high-priority ACMG4_5),
outputting bgzipped and tabix-indexed VCF files with the harmonized naming scheme:
  ACMG{score}.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz
"""

import argparse
import logging
import os
import sys
import pandas as pd
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

ACMG_CLASSES = ["1", "2", "3", "4", "5", "NA", "4_5"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stratify a cohort CNV VCF into ACMG pathogenicity-class VCFs."
    )
    parser.add_argument("--input-vcf", required=True, help="Path to input frequency-filtered cohort VCF (.vcf or .vcf.gz)")
    parser.add_argument("--annotsv-tsv", required=True, help="Path to AnnotSV output TSV file")
    parser.add_argument(
        "--output-pattern",
        default="ACMG{score}.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf.gz",
        help="Output filename pattern with {score} placeholder",
    )
    parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
    return parser.parse_args()


def standardize_chrom(chrom_val):
    c = str(chrom_val).strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def get_acmg_group(val):
    if pd.isna(val):
        return "NA"
    s = str(val).strip()
    if s == "" or s.lower() in ["na", "nan", "<na>"]:
        return "NA"
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
        return str(f)
    except ValueError:
        return s


def load_annotsv_mapping(annotsv_path):
    """
    Returns a dict mapping (std_chrom, pos, end) -> acmg_class (str).
    """
    if not os.path.exists(annotsv_path) or os.path.getsize(annotsv_path) == 0:
        logging.warning(f"AnnotSV file {annotsv_path} is empty or missing.")
        return {}

    try:
        df = pd.read_csv(annotsv_path, sep="\t", keep_default_na=False, low_memory=False)
    except Exception as e:
        logging.error(f"Failed to read AnnotSV TSV: {e}")
        return {}

    if df.empty:
        return {}

    # Prefer 'full' annotation rows representing whole SVs
    if "Annotation_mode" in df.columns:
        df_full = df[df["Annotation_mode"] == "full"]
        if df_full.empty:
            df_full = df
    else:
        df_full = df

    mapping = {}
    for _, row in df_full.iterrows():
        chrom = standardize_chrom(row.get("SV_chrom", ""))
        try:
            pos = int(row.get("SV_start", 0))
            end = int(row.get("SV_end", 0))
        except (ValueError, TypeError):
            continue
        acmg = get_acmg_group(row.get("ACMG_class", "NA"))
        mapping[(chrom, pos, end)] = acmg

    logging.info(f"Loaded {len(mapping)} SV ACMG mappings from AnnotSV.")
    return mapping


def get_record_end(record):
    if "END" in record.info:
        val = record.info["END"]
        if isinstance(val, (tuple, list)):
            return int(val[0])
        return int(val)
    return record.stop


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    acmg_map = load_annotsv_mapping(args.annotsv_tsv)

    try:
        in_vcf = pysam.VariantFile(args.input_vcf)
    except Exception as e:
        logging.error(f"Cannot open input VCF {args.input_vcf}: {e}")
        sys.exit(1)

    header = in_vcf.header
    if "ACMG_CLASS" not in header.info:
        header.info.add("ACMG_CLASS", "1", "String", "AnnotSV ACMG Pathogenicity Class")

    # Initialize record collectors for each score class
    records_by_class = {c: [] for c in ACMG_CLASSES}
    all_annotated_records = []

    for record in in_vcf:
        c_std = standardize_chrom(record.chrom)
        pos = record.pos
        end = get_record_end(record)

        # Look up ACMG class with exact match or fallbacks
        acmg = acmg_map.get((c_std, pos, end))
        if acmg is None:
            # Fallback check matching on chrom & pos
            for (c, p, e), val in acmg_map.items():
                if c == c_std and abs(p - pos) <= 100:
                    acmg = val
                    break
        if acmg is None:
            acmg = "NA"

        try:
            record.info["ACMG_CLASS"] = acmg
        except Exception:
            pass
        all_annotated_records.append(record)

        if acmg in records_by_class:
            records_by_class[acmg].append(record)
        else:
            records_by_class["NA"].append(record)

        # High priority Class 4 & 5 combination
        if acmg in ["4", "5"]:
            records_by_class["4_5"].append(record)

    # 1. Output unified Variant Annotated CNV Callset (all records with ACMG annotation)
    annotated_raw = os.path.join(args.outdir, "annotated.frequency-filtered.cohort.intra-sample.normalised.quality-filtered.cnv.vcf")
    annotated_gz = annotated_raw + ".gz"
    with pysam.VariantFile(annotated_raw, "w", header=header) as out_f:
        for rec in all_annotated_records:
            out_f.write(rec)
    pysam.tabix_compress(annotated_raw, annotated_gz, force=True)
    try:
        pysam.tabix_index(annotated_gz, preset="vcf", force=True)
    except Exception as e:
        logging.warning(f"tabix index warning for {annotated_gz}: {e}")
    if os.path.exists(annotated_raw):
        os.remove(annotated_raw)
    logging.info(f"Written {len(all_annotated_records)} annotated calls -> {annotated_gz}")

    # 2. Output 7 Stratified ACMG CNV Callsets
    for score in ACMG_CLASSES:
        out_name = args.output_pattern.format(score=score)
        if out_name.endswith(".gz"):
            raw_vcf_name = out_name[:-3]
            final_gz_name = out_name
        else:
            raw_vcf_name = out_name
            final_gz_name = out_name + ".gz"

        raw_path = os.path.join(args.outdir, raw_vcf_name)
        final_gz_path = os.path.join(args.outdir, final_gz_name)

        # Write uncompressed VCF
        with pysam.VariantFile(raw_path, "w", header=header) as out_f:
            for rec in records_by_class[score]:
                out_f.write(rec)

        # Bgzip compress and tabix index
        pysam.tabix_compress(raw_path, final_gz_path, force=True)
        try:
            pysam.tabix_index(final_gz_path, preset="vcf", force=True)
        except Exception as e:
            logging.warning(f"tabix index warning for {final_gz_path}: {e}")

        # Clean up temporary raw uncompressed VCF
        if os.path.exists(raw_path):
            os.remove(raw_path)

        logging.info(
            f"Written {len(records_by_class[score])} calls to ACMG class {score} -> {final_gz_path}"
        )


if __name__ == "__main__":
    main()
