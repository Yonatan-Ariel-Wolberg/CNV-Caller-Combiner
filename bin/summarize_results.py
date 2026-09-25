#!/usr/bin/env python3
"""
summarize_results.py

Beginner-friendly results aggregator for CNV-Caller-Combiner.
Scans the pipeline output directory, tallies ACMG clinical classifications,
identifies candidate pathogenic loci, and generates an executive summary
report in both Markdown and terminal console formats.
"""

import argparse
import glob
import os
import sys
from pathlib import Path


def summarize_acmg_files(outdir: Path) -> dict:
    summary = {
        "ACMG_5": 0,
        "ACMG_4": 0,
        "ACMG_4_5": 0,
        "ACMG_3": 0,
        "ACMG_2": 0,
        "ACMG_1": 0,
        "ACMG_NA": 0,
        "total": 0,
        "pathogenic_genes": set(),
        "sample_counts": {},
        "files_found": {}
    }

    # Search for stratified TSVs/CSVs
    tsvs = list(outdir.rglob("*ACMG*.tsv")) + list(outdir.rglob("*ACMG*.csv"))
    for tsv in tsvs:
        stem = tsv.stem
        # Identify tier
        tier = None
        for candidate in ["ACMG_4_5", "ACMG_5", "ACMG_4", "ACMG_3", "ACMG_2", "ACMG_1", "ACMG_NA"]:
            if candidate in stem:
                tier = candidate
                break
        if not tier:
            continue

        summary["files_found"][tier] = str(tsv)
        try:
            with open(tsv, errors="ignore") as f:
                header_line = f.readline()
                header = header_line.strip().split("\t" if "\t" in header_line else ",")
                gene_col = None
                sample_col = None
                for i, col in enumerate(header):
                    col_lower = col.lower()
                    if col_lower in ["genes", "gene", "gene_name", "symbol"]:
                        gene_col = i
                    if col_lower in ["sample", "sampleid", "samples"]:
                        sample_col = i

                count = 0
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.strip().split("\t" if "\t" in line else ",")
                    count += 1
                    if tier in ["ACMG_5", "ACMG_4", "ACMG_4_5"]:
                        if gene_col is not None and len(parts) > gene_col:
                            genes = [g.strip() for g in parts[gene_col].split(";") if g.strip()]
                            summary["pathogenic_genes"].update(genes)
                        if sample_col is not None and len(parts) > sample_col:
                            s = parts[sample_col]
                            summary["sample_counts"][s] = summary["sample_counts"].get(s, 0) + 1

                summary[tier] = count
        except Exception:
            pass

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate and summarize CNV-Caller-Combiner output results."
    )
    parser.add_argument(
        "--results-dir", "-r",
        required=True,
        help="Directory containing the pipeline results (e.g. /path/to/outdir/balanced)."
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Path to save the Markdown summary report (default: <results-dir>/CLINICAL_SUMMARY_REPORT.md)."
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    if not results_dir.is_dir():
        sys.stderr.write(f"❌ Error: Results directory does not exist: {results_dir}\n")
        sys.exit(1)

    print(f"📊 Analyzing pipeline results in: {results_dir} ...")
    summary = summarize_acmg_files(results_dir)

    out_file = Path(args.output) if args.output else (results_dir / "CLINICAL_SUMMARY_REPORT.md")

    report_lines = [
        "# 🧬 CNV-Caller-Combiner: Executive Summary Report",
        "",
        f"**Run Directory:** `{results_dir}`  ",
        f"**Date Generated:** `{os.popen('date').read().strip()}`  ",
        "",
        "---",
        "",
        "## 1. Clinical ACMG Classification Breakdown",
        "",
        "The cohort consensus variants were filtered and stratified into standard clinical tiers according to ACMG/ClinGen guidelines:",
        "",
        "| ACMG Tier | Clinical Significance | Variant Count | Key File Location |",
        "| :--- | :--- | :--- | :--- |",
        f"| **ACMG Class 5** | 🔴 **Pathogenic** | **{summary['ACMG_5']}** | `{summary['files_found'].get('ACMG_5', 'Pending / Under Process')}` |",
        f"| **ACMG Class 4** | 🟠 **Likely Pathogenic** | **{summary['ACMG_4']}** | `{summary['files_found'].get('ACMG_4', 'Pending / Under Process')}` |",
        f"| **ACMG Classes 4 & 5** | 🚨 **High-Confidence Diagnostic** | **{summary['ACMG_4_5']}** | `{summary['files_found'].get('ACMG_4_5', 'Pending / Under Process')}` |",
        f"| **ACMG Class 3** | 🟡 **Variant of Uncertain Significance (VUS)** | **{summary['ACMG_3']}** | `{summary['files_found'].get('ACMG_3', 'Pending / Under Process')}` |",
        f"| **ACMG Class 2** | 🟢 Likely Benign | {summary['ACMG_2']} | `{summary['files_found'].get('ACMG_2', 'Pending / Under Process')}` |",
        f"| **ACMG Class 1** | ⚪ Benign | {summary['ACMG_1']} | `{summary['files_found'].get('ACMG_1', 'Pending / Under Process')}` |",
        f"| **ACMG Class NA** | ❔ Unclassified | {summary['ACMG_NA']} | `{summary['files_found'].get('ACMG_NA', 'Pending / Under Process')}` |",
        "",
        "---",
        "",
        "## 2. Top Pathogenic Candidate Genes Affected",
        "",
    ]

    if summary["pathogenic_genes"]:
        top_genes = sorted(list(summary["pathogenic_genes"]))[:25]
        report_lines.append(f"Identified **{len(summary['pathogenic_genes'])}** unique genes in pathogenic/likely-pathogenic loci:")
        report_lines.append("")
        report_lines.append(", ".join([f"`{g}`" for g in top_genes]))
        if len(summary["pathogenic_genes"]) > 25:
            report_lines.append(f"\n*(and {len(summary['pathogenic_genes']) - 25} more...)*")
        report_lines.append("")
    else:
        report_lines.append("*(No pathogenic genes identified or annotation step in progress)*\n")

    report_lines.extend([
        "---",
        "",
        "## 3. How to View Your Results",
        "",
        "1. **In Microsoft Excel or LibreOffice**:",
        f"   - Open the primary high-confidence file: `{summary['files_found'].get('ACMG_4_5', 'consensus_evaluate/*ACMG_4_5.tsv')}`",
        "2. **In IGV (Integrative Genomics Viewer)**:",
        f"   - Load the refined consensus VCF: `{results_dir}/consensus_refined/*.vcf.gz`",
        f"   - Load your patient BAM: `{results_dir}/../converted_bams/*.bam`",
        "3. **In Web Browser**:",
        f"   - Open the interactive MultiQC report: `{results_dir}/multiqc_report.html`",
        "",
        "---",
        "*Report generated automatically by CNV-Caller-Combiner.*"
    ])

    report_text = "\n".join(report_lines)
    with open(out_file, "w") as f:
        f.write(report_text + "\n")

    print(f"✅ Summary report successfully written to: {out_file}\n")
    print("=" * 60)
    print(f"  🔴 Pathogenic Variants (Class 5):        {summary['ACMG_5']}")
    print(f"  🟠 Likely Pathogenic (Class 4):          {summary['ACMG_4']}")
    print(f"  🚨 Combined Pathogenic Tier (Class 4/5): {summary['ACMG_4_5']}")
    print(f"  🟡 Variants of Uncertain Sig (Class 3):  {summary['ACMG_3']}")
    print("=" * 60)


if __name__ == "__main__":
    main()
