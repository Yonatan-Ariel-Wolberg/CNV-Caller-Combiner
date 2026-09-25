#!/usr/bin/env python3
"""
Convert Mosdepth regions coverage (bed.gz) into caller-ready inputs:
  - xhmm:     Consolidated multi-sample DATA.RD.txt matrix for XHMM
  - canoes:   Per-chromosome read count matrices & sample lists for CANOES
  - clamms:   Per-sample coverage BED for CLAMMS
  - cnvkit:   Per-sample .targetcoverage.cnn for CNVkit

Usage:
  # Generate XHMM DATA.RD.txt from Mosdepth region beds:
  convert_mosdepth_to_caller_inputs.py xhmm --beds sample1.regions.bed.gz sample2.regions.bed.gz --samples S1 S2 -o DATA.RD.txt

  # Generate CANOES read count matrices split by chromosome:
  convert_mosdepth_to_caller_inputs.py canoes --beds sample1.regions.bed.gz sample2.regions.bed.gz --samples S1 S2 -o canoes_counts_dir/

  # Generate CLAMMS coverage BED for single sample:
  convert_mosdepth_to_caller_inputs.py clamms --bed sample1.regions.bed.gz -o sample1.coverage.bed

  # Generate CNVkit target coverage .cnn for single sample:
  convert_mosdepth_to_caller_inputs.py cnvkit --bed sample1.regions.bed.gz -o sample1.targetcoverage.cnn --target-bed probes.bed
"""

import argparse
import gzip
import math
import os
import sys
from collections import defaultdict


def open_file(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_regions_bed(path):
    """
    Yields (chrom, start, end, name, depth) from a Mosdepth regions.bed.gz file.
    Handles 4-column (chr, start, end, depth) and 5-column (chr, start, end, name, depth) formats.
    """
    with open_file(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue

            if len(cols) == 4:
                name = '-'
                try:
                    depth = float(cols[3])
                except ValueError:
                    continue
            else:
                name = cols[3]
                try:
                    depth = float(cols[4])
                except ValueError:
                    try:
                        depth = float(cols[3])
                    except ValueError:
                        continue

            yield chrom, start, end, name, depth


# =============================================================================
# CALLER FORMAT CONVERTERS
# =============================================================================

def convert_to_xhmm(bed_files, sample_ids, output_file):
    """
    Produces XHMM DATA.RD.txt matrix.
    Header: <tab>sample1\tsample2\t...
    Rows: chr:start+1-end\tcov1\tcov2\t...
    """
    if len(bed_files) != len(sample_ids):
        raise ValueError(f"Number of bed files ({len(bed_files)}) does not match sample IDs ({len(sample_ids)})")

    print(f"[XHMM] Processing {len(bed_files)} samples into {output_file}...", file=sys.stderr)

    # 1. Read first file to get canonical targets
    targets = []
    first_bed = bed_files[0]
    for chrom, start, end, _, _ in parse_regions_bed(first_bed):
        # GATK 1-based coordinates
        target_name = f"{chrom}:{start+1}-{end}"
        targets.append(target_name)

    n_targets = len(targets)
    print(f"[XHMM] Found {n_targets} target intervals.", file=sys.stderr)

    # 2. Extract coverage columns for each sample
    sample_coverages = []
    for i, (bed, s_id) in enumerate(zip(bed_files, sample_ids)):
        covs = []
        for _, _, _, _, depth in parse_regions_bed(bed):
            covs.append(depth)
        if len(covs) != n_targets:
            raise ValueError(f"Sample {s_id} ({bed}) has {len(covs)} targets, expected {n_targets}")
        sample_coverages.append(covs)

    # 3. Write DATA.RD.txt
    os.makedirs(os.path.dirname(os.path.abspath(output_file)) or '.', exist_ok=True)
    with open(output_file, 'w') as out_fh:
        # Header starts with a tab or empty field
        out_fh.write("\t" + "\t".join(sample_ids) + "\n")
        for t_idx, target_name in enumerate(targets):
            row_vals = [f"{sample_coverages[s_idx][t_idx]:.2f}" for s_idx in range(len(sample_ids))]
            out_fh.write(f"{target_name}\t" + "\t".join(row_vals) + "\n")

    print(f"[XHMM] Written matrix with {n_targets} targets across {len(sample_ids)} samples to {output_file}", file=sys.stderr)


def convert_to_canoes(bed_files, sample_ids, output_dir, read_length=100):
    """
    Produces CANOES per-chromosome read count matrices and sample lists.
    File: ${output_dir}/${chr}_canoes_reads.txt
    Header: (none)
    Columns: chrom\tstart\tend\tcount_s1\tcount_s2\t...
    File: ${output_dir}/${chr}_sample_list
    """
    if len(bed_files) != len(sample_ids):
        raise ValueError(f"Number of bed files ({len(bed_files)}) does not match sample IDs ({len(sample_ids)})")

    os.makedirs(output_dir, exist_ok=True)
    print(f"[CANOES] Processing {len(bed_files)} samples into {output_dir}...", file=sys.stderr)

    # Parse targets grouped by chromosome
    chr_targets = defaultdict(list)
    for chrom, start, end, _, _ in parse_regions_bed(bed_files[0]):
        chr_targets[chrom].append((start, end))

    # Parse coverage per sample grouped by chromosome
    chr_sample_counts = defaultdict(lambda: [[] for _ in range(len(sample_ids))])

    for s_idx, (bed, s_id) in enumerate(zip(bed_files, sample_ids)):
        for chrom, start, end, _, depth in parse_regions_bed(bed):
            # Estimated read count = round(depth * length / read_length)
            length = max(1, end - start)
            count = int(round((depth * length) / float(read_length)))
            chr_sample_counts[chrom][s_idx].append(count)

    # Write per-chromosome files
    for chrom, targets in chr_targets.items():
        count_file = os.path.join(output_dir, f"{chrom}_canoes_reads.txt")
        slist_file = os.path.join(output_dir, f"{chrom}_sample_list")

        # Write sample list
        with open(slist_file, 'w') as sl_fh:
            for s_id in sample_ids:
                sl_fh.write(f"{s_id}\n")

        # Write counts
        counts = chr_sample_counts[chrom]
        with open(count_file, 'w') as out_fh:
            for t_idx, (start, end) in enumerate(targets):
                row_counts = [str(counts[s_idx][t_idx]) for s_idx in range(len(sample_ids))]
                out_fh.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(row_counts) + "\n")

    print(f"[CANOES] Successfully wrote per-chromosome matrices for {len(chr_targets)} chromosomes.", file=sys.stderr)


def convert_to_clamms(bed_file, output_file):
    """
    Produces CLAMMS coverage BED file.
    Columns: chrom\tstart\tend\tmean_depth
    """
    os.makedirs(os.path.dirname(os.path.abspath(output_file)) or '.', exist_ok=True)
    with open(output_file, 'w') as out_fh:
        for chrom, start, end, _, depth in parse_regions_bed(bed_file):
            out_fh.write(f"{chrom}\t{start}\t{end}\t{depth:.6g}\n")


def convert_to_cnvkit(bed_file, output_file, target_bed=None):
    """
    Produces CNVkit .targetcoverage.cnn file.
    Header: chromosome\tstart\tend\tgene\tdepth\tlog2
    """
    # Load gene names if target_bed provided
    gene_map = {}
    if target_bed and os.path.exists(target_bed):
        with open_file(target_bed, 'rt') as gfh:
            for line in gfh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    key = (parts[0], int(parts[1]), int(parts[2]))
                    gene_map[key] = parts[3]

    os.makedirs(os.path.dirname(os.path.abspath(output_file)) or '.', exist_ok=True)
    with open(output_file, 'w') as out_fh:
        out_fh.write("chromosome\tstart\tend\tgene\tdepth\tlog2\n")
        for chrom, start, end, name, depth in parse_regions_bed(bed_file):
            gene = gene_map.get((chrom, start, end), name if name != '.' else '-')
            if depth > 0.0:
                log2_val = math.log2(depth)
            else:
                log2_val = -20.0
            out_fh.write(f"{chrom}\t{start}\t{end}\t{gene}\t{depth:.6g}\t{log2_val:.6g}\n")


# =============================================================================
# MAIN CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Universal Mosdepth Coverage Converter for CNV Callers")
    subparsers = parser.add_subparsers(dest='caller', required=True)

    # XHMM subparser
    xhmm_parser = subparsers.add_parser('xhmm', help="Generate XHMM DATA.RD.txt")
    xhmm_parser.add_argument('--beds', nargs='+', required=True, help="List of Mosdepth regions.bed.gz files")
    xhmm_parser.add_argument('--samples', nargs='+', required=True, help="List of corresponding sample IDs")
    xhmm_parser.add_argument('-o', '--output', required=True, help="Output DATA.RD.txt path")

    # CANOES subparser
    canoes_parser = subparsers.add_parser('canoes', help="Generate CANOES per-chromosome count tables")
    canoes_parser.add_argument('--beds', nargs='+', required=True, help="List of Mosdepth regions.bed.gz files")
    canoes_parser.add_argument('--samples', nargs='+', required=True, help="List of corresponding sample IDs")
    canoes_parser.add_argument('--read-length', type=int, default=100, help="Mean read length [100]")
    canoes_parser.add_argument('-o', '--output-dir', required=True, help="Output directory for per-chr tables")

    # CLAMMS subparser
    clamms_parser = subparsers.add_parser('clamms', help="Generate CLAMMS .coverage.bed")
    clamms_parser.add_argument('--bed', required=True, help="Single sample Mosdepth regions.bed.gz")
    clamms_parser.add_argument('-o', '--output', required=True, help="Output coverage BED path")

    # CNVkit subparser
    cnvkit_parser = subparsers.add_parser('cnvkit', help="Generate CNVkit .targetcoverage.cnn")
    cnvkit_parser.add_argument('--bed', required=True, help="Single sample Mosdepth regions.bed.gz")
    cnvkit_parser.add_argument('--target-bed', default=None, help="Optional original target BED with gene annotations")
    cnvkit_parser.add_argument('-o', '--output', required=True, help="Output .targetcoverage.cnn path")

    args = parser.parse_args()

    if args.caller == 'xhmm':
        convert_to_xhmm(args.beds, args.samples, args.output)
    elif args.caller == 'canoes':
        convert_to_canoes(args.beds, args.samples, args.output_dir, args.read_length)
    elif args.caller == 'clamms':
        convert_to_clamms(args.bed, args.output)
    elif args.caller == 'cnvkit':
        convert_to_cnvkit(args.bed, args.output, args.target_bed)


if __name__ == '__main__':
    main()
