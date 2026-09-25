#!/usr/bin/env python3
"""
compile_quality_summary.py

Compiles multi-stage CNV callset quality metrics summary reports (*_quality_summary.txt)
into unified Markdown and TSV comparison tables.
"""

import argparse
import glob
import os
import sys

def parse_summary_file(filepath):
    stage = os.path.basename(filepath).replace('_quality_summary.txt', '')
    data = {
        'Stage': stage, 'Calls': '0', 'DEL': '0', 'DUP': '0', 'DEL/DUP': 'N/A',
        'Median_Size_bp': '0', 'Median_Size_DEL_bp': '0', 'Median_Size_DUP_bp': '0',
        'Total_Mb': '0', 'DEL_Mb': '0', 'DUP_Mb': '0',
        'Median_Probes': '0', 'Median_Probes_DEL': '0', 'Median_Probes_DUP': '0',
        'Mean_Calls_Sample': '0', 'Mean_Calls_DEL_Sample': '0', 'Mean_Calls_DUP_Sample': '0'
    }
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Total CNV Calls:'):
                    data['Calls'] = line.split(':', 1)[1].strip()
                elif line.startswith('Deletions (DEL):'):
                    data['DEL'] = line.split(':', 1)[1].split('(')[0].strip()
                elif line.startswith('Duplications (DUP):'):
                    data['DUP'] = line.split(':', 1)[1].split('(')[0].strip()
                elif line.startswith('DEL / DUP Ratio:'):
                    data['DEL/DUP'] = line.split(':', 1)[1].strip()
                elif line.startswith('Median Length (DEL):'):
                    data['Median_Size_DEL_bp'] = line.split(':', 1)[1].replace('bp','').strip()
                elif line.startswith('Median Length (DUP):'):
                    data['Median_Size_DUP_bp'] = line.split(':', 1)[1].replace('bp','').strip()
                elif line.startswith('Median Length:'):
                    data['Median_Size_bp'] = line.split(':', 1)[1].replace('bp','').strip()
                elif line.startswith('Total Genomic Footprint (DEL):'):
                    data['DEL_Mb'] = line.split(':', 1)[1].replace('Mb','').strip()
                elif line.startswith('Total Genomic Footprint (DUP):'):
                    data['DUP_Mb'] = line.split(':', 1)[1].replace('Mb','').strip()
                elif line.startswith('Total Genomic Footprint:'):
                    data['Total_Mb'] = line.split(':', 1)[1].replace('Mb','').strip()
                elif line.startswith('Median Probes per CNV (DEL):'):
                    data['Median_Probes_DEL'] = line.split(':', 1)[1].strip()
                elif line.startswith('Median Probes per CNV (DUP):'):
                    data['Median_Probes_DUP'] = line.split(':', 1)[1].strip()
                elif line.startswith('Median Probes per CNV:'):
                    data['Median_Probes'] = line.split(':', 1)[1].strip()
                elif line.startswith('Mean Calls per Sample (DEL):'):
                    data['Mean_Calls_DEL_Sample'] = line.split(':', 1)[1].strip()
                elif line.startswith('Mean Calls per Sample (DUP):'):
                    data['Mean_Calls_DUP_Sample'] = line.split(':', 1)[1].strip()
                elif line.startswith('Mean Calls per Sample:'):
                    data['Mean_Calls_Sample'] = line.split(':', 1)[1].strip()
    except Exception as e:
        sys.stderr.write(f"Warning: could not parse {filepath}: {e}\n")
    return data

def main():
    parser = argparse.ArgumentParser(description="Compile multi-stage CNV quality assessment summaries.")
    parser.add_argument("--pattern", default="*_quality_summary.txt", help="Glob pattern for quality summary files.")
    parser.add_argument("--output-md", default="all_stages_quality_metrics_summary.md", help="Output markdown table path.")
    parser.add_argument("--output-tsv", default="all_stages_quality_metrics_summary.tsv", help="Output TSV table path.")
    args = parser.parse_args()

    reports = sorted(glob.glob(args.pattern))
    rows = [parse_summary_file(r) for r in reports]

    with open(args.output_md, 'w', encoding='utf-8') as f:
        f.write('# Multi-Stage CNV Callset Quality Metrics Assessment\n\n## 1. Overall Callset Metrics & Variant Breakdown\n\n')
        f.write('| Stage | Total Calls | DEL | DUP | DEL/DUP Ratio | Median Size (bp) | Total Footprint (Mb) | Median Probes | Calls / Sample |\n')
        f.write('| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n')
        for r in rows:
            f.write(f"| **{r['Stage']}** | {r['Calls']} | {r['DEL']} | {r['DUP']} | {r['DEL/DUP']} | {r['Median_Size_bp']} | {r['Total_Mb']} | {r['Median_Probes']} | {r['Mean_Calls_Sample']} |\n")

        f.write('\n## 2. 3-Way Comparative Evaluation: Genomic Size & Footprint (Total vs DEL vs DUP)\n\n')
        f.write('| Stage | Total Median (bp) | DEL Median (bp) | DUP Median (bp) | Total Footprint (Mb) | DEL Footprint (Mb) | DUP Footprint (Mb) |\n')
        f.write('| :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n')
        for r in rows:
            f.write(f"| **{r['Stage']}** | {r['Median_Size_bp']} | {r['Median_Size_DEL_bp']} | {r['Median_Size_DUP_bp']} | {r['Total_Mb']} | {r['DEL_Mb']} | {r['DUP_Mb']} |\n")

        f.write('\n## 3. 3-Way Comparative Evaluation: Probe Support & Sample Burden (Total vs DEL vs DUP)\n\n')
        f.write('| Stage | Total Median Probes | DEL Median Probes | DUP Median Probes | Total Calls / Sample | DEL Calls / Sample | DUP Calls / Sample |\n')
        f.write('| :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n')
        for r in rows:
            f.write(f"| **{r['Stage']}** | {r['Median_Probes']} | {r['Median_Probes_DEL']} | {r['Median_Probes_DUP']} | {r['Mean_Calls_Sample']} | {r['Mean_Calls_DEL_Sample']} | {r['Mean_Calls_DUP_Sample']} |\n")

    with open(args.output_tsv, 'w', encoding='utf-8') as f:
        f.write('Stage\tTotal_Calls\tDEL\tDUP\tDEL_DUP_Ratio\tMedian_Size_bp\tMedian_Size_DEL_bp\tMedian_Size_DUP_bp\tTotal_Mb\tDEL_Mb\tDUP_Mb\tMedian_Probes\tMedian_Probes_DEL\tMedian_Probes_DUP\tMean_Calls_Per_Sample\tMean_Calls_DEL_Per_Sample\tMean_Calls_DUP_Per_Sample\n')
        for r in rows:
            f.write(f"{r['Stage']}\t{r['Calls']}\t{r['DEL']}\t{r['DUP']}\t{r['DEL/DUP']}\t{r['Median_Size_bp']}\t{r['Median_Size_DEL_bp']}\t{r['Median_Size_DUP_bp']}\t{r['Total_Mb']}\t{r['DEL_Mb']}\t{r['DUP_Mb']}\t{r['Median_Probes']}\t{r['Median_Probes_DEL']}\t{r['Median_Probes_DUP']}\t{r['Mean_Calls_Sample']}\t{r['Mean_Calls_DEL_Sample']}\t{r['Mean_Calls_DUP_Sample']}\n")

if __name__ == '__main__':
    main()
