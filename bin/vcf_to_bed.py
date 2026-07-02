#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: vcf_to_bed.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
"""
Convert a single-sample CNV VCF file to a 5-column BED file.

Output columns (tab-separated, 0-based half-open coordinates):
    CHR  START  STOP  CNV_TYPE  SAMPLE_ID

CNV_TYPE is derived from the ALT field (e.g. <DEL> → DEL, <DUP> → DUP).
If the ALT field does not resolve to DEL or DUP the value CNV is used.
STOP is taken from the INFO/END field when present; otherwise START+1 is used.
SAMPLE_ID is read from the single-sample column header in the VCF.

Usage:
    vcf_to_bed.py --vcf sample.sorted.vcf.gz --output sample.bed
    vcf_to_bed.py --vcf sample.vcf --output sample.bed [--sample_id OVERRIDE]
"""

import argparse
import gzip
import os
import sys


def _open_vcf(path):
    """Return a file handle that reads the VCF regardless of compression."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def vcf_to_bed(vcf_path, output_path, sample_id_override=None):
    """Parse *vcf_path* and write a 5-column BED to *output_path*."""
    sample_id = sample_id_override
    records = []

    with _open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Skip meta-information lines
            if line.startswith('##'):
                continue

            # Column-header line – extract sample name from column 10
            if line.startswith('#CHROM'):
                if sample_id is None:
                    cols = line.split('\t')
                    if len(cols) > 9:
                        sample_id = cols[9]
                continue

            if not line:
                continue

            cols = line.split('\t')
            if len(cols) < 5:
                continue  # skip malformed lines

            chrom = cols[0]
            pos   = int(cols[1])          # 1-based
            start = pos - 1               # convert to 0-based
            alt   = cols[4].strip('<>')   # e.g. <DEL> → DEL
            info  = cols[7] if len(cols) > 7 else ''

            # Derive CNV type from ALT field
            alt_upper = alt.upper()
            if alt_upper in ('DEL', 'DUP'):
                cnv_type = alt_upper
            else:
                cnv_type = 'CNV'

            # Extract END coordinate from INFO field
            end = None
            for token in info.split(';'):
                if token.startswith('END='):
                    try:
                        end = int(token.split('=', 1)[1])
                    except ValueError:
                        pass
                    break
            if end is None:
                end = start + 1  # fallback when END is absent

            records.append((chrom, start, end, cnv_type))

    if sample_id is None:
        sample_id = os.path.basename(vcf_path).split('.')[0]

    with open(output_path, 'w') as out:
        for chrom, start, end, cnv_type in records:
            out.write(f'{chrom}\t{start}\t{end}\t{cnv_type}\t{sample_id}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Convert a single-sample CNV VCF to a 5-column BED file.')
    parser.add_argument('--vcf',       required=True,
                        help='Input VCF or VCF.gz file path')
    parser.add_argument('--output',    required=True,
                        help='Output BED file path')
    parser.add_argument('--sample_id', default=None,
                        help='Override sample ID (default: read from VCF header)')
    args = parser.parse_args()

    vcf_to_bed(args.vcf, args.output, sample_id_override=args.sample_id)


if __name__ == '__main__':
    main()
