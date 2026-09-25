#!/usr/bin/env python3
"""
make_samplesheet.py

Automated, beginner-friendly samplesheet generator for CNV-Caller-Combiner.
Scans directories for BAM or CRAM alignment files, validates index presence,
sanitizes sample identifiers, and outputs a formatted TSV samplesheet.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path


def clean_sample_id(filename: str) -> str:
    """Derive a clean SampleID from a file name."""
    name = Path(filename).name
    name = re.sub(r'\.(bam|cram)(\.(bai|crai))?$', '', name, flags=re.I)
    name = re.sub(r'\.md$', '', name, flags=re.I)
    name = re.sub(r'\.(sorted|aligned|dedup|recal|raw)$', '', name, flags=re.I)
    return name


def find_index(align_path: Path, ext_type: str) -> Path | None:
    """Find the index file corresponding to a BAM or CRAM file."""
    if ext_type == "bam":
        candidates = [
            align_path.with_suffix(align_path.suffix + ".bai"),
            align_path.with_suffix(".bai"),
        ]
    else:
        candidates = [
            align_path.with_suffix(align_path.suffix + ".crai"),
            align_path.with_suffix(".crai"),
        ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Auto-generate a validated samplesheet from a folder of BAMs or CRAMs."
    )
    parser.add_argument(
        "--input-dir", "-i",
        required=True,
        help="Directory containing BAM or CRAM alignment files (scanned recursively)."
    )
    parser.add_argument(
        "--output", "-o",
        default="samplesheet.tsv",
        help="Output TSV samplesheet path (default: samplesheet.tsv)."
    )
    parser.add_argument(
        "--type", "-t",
        choices=["auto", "bam", "cram"],
        default="auto",
        help="Alignment file type to search for (default: auto-detect)."
    )
    parser.add_argument(
        "--max-samples", "-m",
        type=int,
        default=None,
        help="Optional limit on the number of samples to include (useful for testing)."
    )

    args = parser.parse_args()
    input_dir = Path(args.input_dir).resolve()

    if not input_dir.is_dir():
        sys.stderr.write(f"❌ Error: Input directory does not exist: {input_dir}\n")
        sys.exit(1)

    print(f"🔍 Scanning directory: {input_dir} ...")

    bams = sorted(input_dir.rglob("*.bam")) if args.type in ["auto", "bam"] else []
    crams = sorted(input_dir.rglob("*.cram")) if args.type in ["auto", "cram"] else []

    if args.type == "auto":
        if crams and not bams:
            chosen_type = "cram"
            files = crams
        elif bams and not crams:
            chosen_type = "bam"
            files = bams
        elif crams and bams:
            if len(bams) >= len(crams):
                chosen_type = "bam"
                files = bams
            else:
                chosen_type = "cram"
                files = crams
        else:
            sys.stderr.write(f"❌ No .bam or .cram files found inside {input_dir}\n")
            sys.exit(1)
    elif args.type == "bam":
        chosen_type = "bam"
        files = bams
    else:
        chosen_type = "cram"
        files = crams

    if not files:
        sys.stderr.write(f"❌ No .{chosen_type} files found in {input_dir}\n")
        sys.exit(1)

    print(f"✅ Found {len(files)} {chosen_type.upper()} file(s). Validating index pairs...")

    rows = []
    missing_indices = []
    seen_ids = set()

    for f in files:
        idx = find_index(f, chosen_type)
        sid = clean_sample_id(f.name)

        if sid in seen_ids:
            sid = f"{sid}_{len(seen_ids)+1}"
        seen_ids.add(sid)

        if not idx:
            missing_indices.append(f)
            continue

        rows.append((sid, str(f.resolve()), str(idx.resolve())))

    if missing_indices:
        print(f"⚠️ Warning: {len(missing_indices)} file(s) missing an index (.{chosen_type}.bai / .{chosen_type}.crai):")
        for m in missing_indices[:5]:
            print(f"   - {m.name}")
        if len(missing_indices) > 5:
            print(f"   ... and {len(missing_indices) - 5} more.")
        print("💡 Hint: Generate missing indices with: samtools index <file>")

    if not rows:
        sys.stderr.write("❌ No valid alignment/index pairs found. Cannot create samplesheet.\n")
        sys.exit(1)

    if args.max_samples and args.max_samples > 0:
        rows = rows[:args.max_samples]
        print(f"ℹ️ Capped at {len(rows)} sample(s) as requested.")

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = "SampleID\tBAM\tBAI\n" if chosen_type == "bam" else "SampleID\tCRAM\tCRAI\n"
    with open(out_path, "w") as out:
        out.write(header)
        for sid, align, idx in rows:
            out.write(f"{sid}\t{align}\t{idx}\n")

    print(f"🎉 Successfully created samplesheet with {len(rows)} samples!")
    print(f"📁 Saved to: {out_path}\n")
    print("Top 3 entries preview:")
    for sid, align, idx in rows[:3]:
        print(f"   • {sid} -> {Path(align).name} (Index: {Path(idx).name})")
    print()

if __name__ == "__main__":
    main()
