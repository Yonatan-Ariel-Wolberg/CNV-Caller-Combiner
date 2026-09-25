#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: vcf_to_bed.py
# Purpose: This pipeline helper script converts a single-sample CNV VCF file to 
#   a standard 5-column BED file format.
# Scope: General utility execution.
# Inputs:
#   - --vcf: Input VCF or VCF.gz file path containing CNV caller output.
#   - --sample_id: Optional sample ID to override the name extracted from VCF header.
# Outputs:
#   - --output: Path to the generated 5-column tab-separated BED file.
# Operational notes:
#   - Coordinates in VCF are 1-based closed, while BED uses 0-based half-open.
#     This script adjusts the start coordinate accordingly (start = POS - 1).
#   - CNV type is extracted and normalized from the ALT field (e.g. <DEL> -> DEL).
#   - END coordinate is parsed from the INFO/END tag. If absent, defaults to POS.
#   - Keeps parameter names and structure synchronized with pipeline usage.
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
    """
    Return a file handle that reads the VCF regardless of compression.

    Parameters:
    - path: File path to VCF (automatically checks for .gz suffix).

    Returns:
    - Text file reader context.
    """
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def vcf_to_bed(vcf_path, output_path, sample_id_override=None):
    """
    Parse *vcf_path* and write a 5-column BED to *output_path*.

    Parameters:
    - vcf_path: Path to the input VCF/VCF.gz file.
    - output_path: Path to write the output BED file.
    - sample_id_override: Optional sample ID override string.
    """
    sample_id = sample_id_override
    records = []

    if not os.path.exists(vcf_path) or os.path.getsize(vcf_path) == 0:
        if sample_id is None:
            sample_id = os.path.basename(vcf_path).split('.')[0]
        with open(output_path, 'w') as out:
            pass
        return

    with _open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Skip VCF meta-information header lines
            if line.startswith('##'):
                continue

            # Column-header line – extract sample name from the 10th column (index 9)
            if line.startswith('#CHROM'):
                if sample_id is None:
                    cols = line.split('\t')
                    if len(cols) > 9:
                        sample_id = cols[9]
                continue

            # Skip empty lines
            if not line:
                continue

            cols = line.split('\t')
            if len(cols) < 5:
                continue  # skip malformed lines

            chrom = cols[0]
            try:
                pos = int(cols[1])          # 1-based coordinate in VCF
            except (ValueError, TypeError):
                continue
            start = max(0, pos - 1)         # convert to 0-based coordinate for BED (invariant: start >= 0)
            alt   = cols[4].strip('<>')   # strip brackets to normalize, e.g. <DEL> → DEL
            info  = cols[7] if len(cols) > 7 else ''
            # Normalize and derive CNV type from ALT allele field
            alt_upper = alt.upper()
            if alt_upper in ('DEL', 'DUP') or alt_upper.startswith('DEL') or alt_upper.startswith('DUP'):
                cnv_type = 'DEL' if 'DEL' in alt_upper else 'DUP'
            elif alt_upper.startswith('CN'):
                try:
                    cn_val = int(alt_upper.replace('CN', ''))
                    if cn_val < 2:
                        cnv_type = 'DEL'
                    elif cn_val > 2:
                        cnv_type = 'DUP'
                    else:
                        cnv_type = 'CNV'
                except:
                    cnv_type = 'CNV'
            else:
                cnv_type = 'CNV'
            # Extract END coordinate from INFO field by parsing key-value pairs
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
            else:
                end = max(start + 1, end)  # invariant: end > start

            records.append((chrom, start, end, cnv_type))

    # Fallback to file basename prefix if sample ID could not be determined
    if sample_id is None:
        sample_id = os.path.basename(vcf_path).split('.')[0]

    # Write the records out to the output path in tab-separated BED format
    with open(output_path, 'w') as out:
        for chrom, start, end, cnv_type in records:
            out.write(f'{chrom}\t{start}\t{end}\t{cnv_type}\t{sample_id}\n')


def main():
    """Parse command line arguments and execute VCF-to-BED conversion."""
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

