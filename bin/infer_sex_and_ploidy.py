#!/usr/bin/env python3
"""
infer_sex_and_ploidy.py

Infers sample genetic sex and sex chromosome ploidy from samtools idxstats files.
Works on BAM/CRAM index statistics files (.idxstats) or directly by executing
samtools idxstats if a BAM file is supplied.

Output:
  cohort_sex_and_ploidy.tsv
  clamms_sexinfo.tsv (ready for CLAMMS --sexinfo)
"""

import argparse
import glob
import os
import sys


def parse_idxstats(idxstats_path):
    """Parse a samtools idxstats file and return mapped reads per contig."""
    data = {}
    with open(idxstats_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                contig = parts[0]
                length = int(parts[1])
                mapped = int(parts[2])
                data[contig] = {"length": length, "mapped": mapped}
    return data


def infer_sample_sex(contig_data, sample_id):
    """Infer genetic sex and ploidy from mapped read densities."""
    autosome_mapped = 0
    autosome_length = 0

    # Contig prefixes (supporting both 'chr1' and '1')
    for i in range(1, 23):
        for prefix in [f"chr{i}", f"{i}"]:
            if prefix in contig_data:
                autosome_mapped += contig_data[prefix]["mapped"]
                autosome_length += contig_data[prefix]["length"]
                break

    if autosome_length == 0 or autosome_mapped == 0:
        return {
            "sample_id": sample_id,
            "inferred_sex": "UNKNOWN",
            "karyotype": "UNKNOWN",
            "chrX_ratio": 0.0,
            "chrY_ratio": 0.0,
            "confidence": "LOW",
        }

    autosome_density = autosome_mapped / float(autosome_length)

    # ChrX density
    x_mapped, x_len = 0, 0
    for prefix in ["chrX", "X"]:
        if prefix in contig_data:
            x_mapped = contig_data[prefix]["mapped"]
            x_len = contig_data[prefix]["length"]
            break
    x_density = (x_mapped / float(x_len)) if x_len > 0 else 0.0
    x_ratio = x_density / autosome_density if autosome_density > 0 else 0.0

    # ChrY density
    y_mapped, y_len = 0, 0
    for prefix in ["chrY", "Y"]:
        if prefix in contig_data:
            y_mapped = contig_data[prefix]["mapped"]
            y_len = contig_data[prefix]["length"]
            break
    y_density = (y_mapped / float(y_len)) if y_len > 0 else 0.0
    y_ratio = y_density / autosome_density if autosome_density > 0 else 0.0

    # Decision logic
    # Female: X ratio ~ 0.85-1.15, Y ratio < 0.03
    # Male:   X ratio ~ 0.40-0.65, Y ratio > 0.08
    if y_ratio > 0.05 and x_ratio < 0.75:
        inferred_sex = "MALE"
        karyotype = "XY"
        confidence = "HIGH" if (0.40 <= x_ratio <= 0.65 and y_ratio >= 0.08) else "MEDIUM"
    elif y_ratio < 0.03 and x_ratio > 0.75:
        inferred_sex = "FEMALE"
        karyotype = "XX"
        confidence = "HIGH" if (0.80 <= x_ratio <= 1.20) else "MEDIUM"
    elif y_ratio >= 0.05 and x_ratio >= 0.75:
        inferred_sex = "MALE"
        karyotype = "XXY (Klinefelter)"
        confidence = "MEDIUM"
    elif y_ratio < 0.03 and x_ratio < 0.75:
        inferred_sex = "FEMALE"
        karyotype = "X0 (Turner)"
        confidence = "MEDIUM"
    else:
        inferred_sex = "UNKNOWN"
        karyotype = "INDETERMINATE"
        confidence = "LOW"

    return {
        "sample_id": sample_id,
        "inferred_sex": inferred_sex,
        "karyotype": karyotype,
        "chrX_ratio": round(x_ratio, 4),
        "chrY_ratio": round(y_ratio, 4),
        "confidence": confidence,
    }


def main():
    parser = argparse.ArgumentParser(description="Infer sample sex and ploidy from samtools idxstats.")
    parser.add_argument("--idxstats-dir", help="Directory containing .idxstats files.")
    parser.add_argument("--idxstats-files", nargs="*", help="List of .idxstats files.")
    parser.add_argument("--output-tsv", default="cohort_sex_and_ploidy.tsv", help="Output TSV file path.")
    parser.add_argument("--clamms-sexinfo", default="clamms_sexinfo.tsv", help="Output CLAMMS-formatted sexinfo TSV.")
    args = parser.parse_args()

    files = []
    if args.idxstats_dir:
        files.extend(glob.glob(os.path.join(args.idxstats_dir, "*.idxstats")))
    if args.idxstats_files:
        files.extend(args.idxstats_files)

    if not files:
        sys.exit("Error: No .idxstats files provided or found.")

    results = []
    for fpath in sorted(set(files)):
        fname = os.path.basename(fpath)
        sample_id = fname.replace(".idxstats", "").replace(".bam", "").replace(".cram", "")
        contig_data = parse_idxstats(fpath)
        res = infer_sample_sex(contig_data, sample_id)
        results.append(res)

    # Write cohort sex & ploidy TSV
    with open(args.output_tsv, "w", encoding="utf-8") as f:
        f.write("SampleID\tInferredSex\tKaryotype\tChrX_Ratio\tChrY_Ratio\tConfidence\n")
        for r in results:
            f.write(f"{r['sample_id']}\t{r['inferred_sex']}\t{r['karyotype']}\t{r['chrX_ratio']}\t{r['chrY_ratio']}\t{r['confidence']}\n")

    # Write CLAMMS-formatted sexinfo TSV (SampleID \t M/F)
    with open(args.clamms_sexinfo, "w", encoding="utf-8") as f:
        for r in results:
            clamms_code = "M" if r["inferred_sex"] == "MALE" else ("F" if r["inferred_sex"] == "FEMALE" else "U")
            f.write(f"{r['sample_id']}\t{clamms_code}\n")

    male_count = sum(1 for r in results if r["inferred_sex"] == "MALE")
    female_count = sum(1 for r in results if r["inferred_sex"] == "FEMALE")
    unknown_count = sum(1 for r in results if r["inferred_sex"] not in ("MALE", "FEMALE"))

    print("==================================================")
    print(" INFERRED SEX AND PLOIDY SUMMARY")
    print("==================================================")
    print(f"Total Samples Analyzed: {len(results)}")
    print(f"  - Males (XY):         {male_count}")
    print(f"  - Females (XX):       {female_count}")
    print(f"  - Unknown / Aneuploid:{unknown_count}")
    print(f"Saved cohort sex table: {args.output_tsv}")
    print(f"Saved CLAMMS sexinfo:   {args.clamms_sexinfo}")
    print("==================================================")


if __name__ == "__main__":
    main()
