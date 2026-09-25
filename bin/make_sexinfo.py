#!/usr/bin/env python3
"""
make_sexinfo.py

Automated sexinfo generator for CNV-Caller-Combiner.
Generates a valid tab-delimited sample sex metadata table (sample_id, sex)
required by CLAMMS and GATK ploidy callers.

Capabilities:
1. Derives sample IDs directly from samplesheet.tsv or alignment folder.
2. Inferred biological sex from chromosome X and Y read counts via 'samtools idxstats' when available.
3. Provides a clean fallback default sex ('F' or 'M') if alignments cannot be scanned.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def infer_sex_from_bam(bam_path: str) -> str | None:
    """Infer chromosomal sex (M or F) from chrX and chrY mapped read counts via samtools idxstats."""
    if not shutil.which("samtools") or not os.path.exists(bam_path):
        return None
    try:
        res = subprocess.run(
            ["samtools", "idxstats", bam_path],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if res.returncode != 0:
            return None

        x_reads = 0
        y_reads = 0
        for line in res.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                rname = parts[0].lower().replace("chr", "")
                mapped = int(parts[2])
                if rname == "x":
                    x_reads += mapped
                elif rname == "y":
                    y_reads += mapped

        if x_reads + y_reads == 0:
            return None

        # Y/X read ratio threshold: male typically has Y/X ratio > 0.05 in WES, female < 0.01
        y_ratio = y_reads / float(x_reads + 1)
        return "M" if y_ratio >= 0.03 else "F"
    except Exception:
        return None


def parse_samplesheet(samplesheet_path: str) -> list[tuple[str, str]]:
    """Extract (sample_id, filepath) pairs from samplesheet TSV."""
    samples = []
    with open(samplesheet_path, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        id_col = 0
        path_col = 1
        found_alignment = False
        for idx, h in enumerate(header):
            h_lower = h.lower()
            if "sample" in h_lower or h_lower == "id":
                id_col = idx
            elif any(x in h_lower for x in ["bai", "crai", "index"]):
                continue
            elif any(x in h_lower for x in ["bam", "cram"]):
                path_col = idx
                found_alignment = True
            elif not found_alignment and any(x in h_lower for x in ["file", "path"]):
                path_col = idx

        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > id_col:
                s_id = parts[id_col].strip()
                s_path = parts[path_col].strip() if len(parts) > path_col else ""
                if s_id:
                    samples.append((s_id, s_path))
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Auto-generate a valid sexinfo table from a samplesheet or alignment list."
    )
    parser.add_argument(
        "--samplesheet", "-s",
        required=True,
        help="Path to samplesheet.tsv."
    )
    parser.add_argument(
        "--output", "-o",
        default="sexinfo.tsv",
        help="Output sexinfo TSV file path (default: sexinfo.tsv)."
    )
    parser.add_argument(
        "--default-sex", "-d",
        default="F",
        choices=["M", "F"],
        help="Default sex assigned when alignments cannot be probed (default: F)."
    )
    parser.add_argument(
        "--auto-infer",
        action="store_true",
        default=True,
        help="Attempt biological sex inference via samtools idxstats (default: True)."
    )
    args = parser.parse_args()

    if not os.path.exists(args.samplesheet):
        sys.exit(f"Error: Samplesheet does not exist: {args.samplesheet}")

    samples = parse_samplesheet(args.samplesheet)
    if not samples:
        sys.exit(f"Error: No sample entries found in {args.samplesheet}")

    out_lines = ["sample_id\tsex\n"]
    inferred_count = 0
    fallback_count = 0

    for s_id, s_path in samples:
        sex = None
        if args.auto_infer and s_path and s_path.endswith(".bam") and os.path.exists(s_path):
            sex = infer_sex_from_bam(s_path)
            if sex:
                inferred_count += 1

        if not sex:
            sex = args.default_sex
            fallback_count += 1

        out_lines.append(f"{s_id}\t{sex}\n")

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(out_lines)

    print(f"Generated sexinfo file for {len(samples)} samples: {out_path}")
    print(f"  - Inferred via coverage (X/Y): {inferred_count}")
    print(f"  - Default/Fallback assigned:   {fallback_count} (Sex: {args.default_sex})")


if __name__ == "__main__":
    main()
