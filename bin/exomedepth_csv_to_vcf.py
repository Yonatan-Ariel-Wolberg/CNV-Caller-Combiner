#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: exomedepth_csv_to_vcf.py
# Purpose: Converts tabular CSV CNV calls from the ExomeDepth caller into
#   standard single-sample VCF files suitable for multi-caller merging (Truvari/SURVIVOR).
#
# Arguments:
#   - --input_file: Path to ExomeDepth output CSV file.
#   - --sample_id: Sample identifier corresponding to the analyzed sample.
#   - --output_dir: Target directory path to write the VCF file.
#   - --fai_file: Path to reference FASTA index (.fai) for contig headers.
#   - --min_bf: (Optional) Minimum Bayes Factor for PASS filter. Default: 10.0.
#
# Output:
#   - <output_dir>/<sample_id>_EXOMEDEPTH_output.vcf
# =============================================================================

import argparse
import os
import csv
import logging
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="Convert ExomeDepth CSV calls to standard VCF format.")
    parser.add_argument("--input_file", required=True, help="Path to input ExomeDepth CSV file.")
    parser.add_argument("--sample_id", required=True, help="Sample ID.")
    parser.add_argument("--output_dir", required=True, help="Output directory to save VCF file.")
    parser.add_argument("--fai_file", required=True, help="Path to FASTA index (.fai) reference file.")
    parser.add_argument("--min_bf", type=float, default=10.0, help="Minimum Bayes Factor threshold for PASS (default: 10.0).")
    return parser.parse_args()


def load_fai_contigs(fai_path):
    """
    Parse FASTA index (.fai) to extract chromosome names and lengths for VCF header.
    """
    contigs = []
    if not os.path.exists(fai_path):
        logging.warning(f"FAI file {fai_path} not found. Skipping contig header lines.")
        return contigs

    with open(fai_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                length = parts[1]
                contigs.append((chrom, length))
    return contigs


def write_vcf_header(vcf_out, contigs, sample_id):
    """
    Write standard VCFv4.2 headers with necessary INFO, FILTER, and FORMAT definitions.
    """
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write(f"##fileDate={datetime.today().strftime('%Y%m%d')}\n")
    vcf_out.write("##source=ExomeDepth\n")

    # Contigs
    for chrom, length in contigs:
        vcf_out.write(f"##contig=<ID={chrom},length={length}>\n")

    # ALT tags
    vcf_out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    vcf_out.write('##ALT=<ID=DUP,Description="Duplication">\n')

    # FILTER tags
    vcf_out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    vcf_out.write('##FILTER=<ID=LowQual,Description="Low confidence Bayes Factor score">\n')

    # INFO tags
    vcf_out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    vcf_out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">\n')
    vcf_out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    vcf_out.write('##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool name">\n')
    vcf_out.write('##INFO=<ID=BF,Number=1,Type=Float,Description="ExomeDepth Bayes Factor">\n')
    vcf_out.write('##INFO=<ID=READS_RATIO,Number=1,Type=Float,Description="Ratio of observed over expected read count">\n')
    vcf_out.write('##INFO=<ID=READS_EXPECTED,Number=1,Type=Float,Description="Expected read count under normal diploid model">\n')
    vcf_out.write('##INFO=<ID=READS_OBSERVED,Number=1,Type=Float,Description="Observed read count in test sample">\n')
    vcf_out.write('##INFO=<ID=NEXONS,Number=1,Type=Integer,Description="Number of target exons spanned by CNV">\n')

    # FORMAT tags
    vcf_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_out.write('##FORMAT=<ID=BF,Number=1,Type=Float,Description="Bayes Factor quality metric">\n')
    vcf_out.write('##FORMAT=<ID=RR,Number=1,Type=Float,Description="Ratio of observed to expected reads">\n')
    vcf_out.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Estimated copy number">\n')

    # Column header
    vcf_out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")


def convert_csv_to_vcf(input_csv, sample_id, output_dir, fai_file, min_bf):
    os.makedirs(output_dir, exist_ok=True)
    out_vcf_path = os.path.join(output_dir, f"{sample_id}_EXOMEDEPTH_output.vcf")

    contigs = load_fai_contigs(fai_file)

    records = []
    if os.path.exists(input_csv) and os.path.getsize(input_csv) > 0:
        with open(input_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Expected fields: chromosome, start, end, type, BF, reads.expected, reads.observed, reads.ratio, nexons
                chrom = row.get("chromosome", "").strip()
                if not chrom:
                    continue

                try:
                    start = int(row["start"])
                    end = int(row["end"])
                    bf = float(row.get("BF", 0.0))
                    reads_ratio = float(row.get("reads.ratio", 1.0))
                    reads_exp = float(row.get("reads.expected", 0.0))
                    reads_obs = float(row.get("reads.observed", 0.0))
                    nexons = int(row.get("nexons", 1))
                except (ValueError, KeyError) as e:
                    logging.warning(f"Skipping malformed ExomeDepth row in {input_csv}: {row} ({e})")
                    continue

                call_type = row.get("type", "").strip().lower()
                if "del" in call_type:
                    svtype = "DEL"
                    alt = "<DEL>"
                    svlen = -(end - start + 1)
                    # Infer genotype and copy number from reads ratio
                    if reads_ratio < 0.15:
                        cn = 0
                        gt = "1/1"
                    else:
                        cn = 1
                        gt = "0/1"
                elif "dup" in call_type:
                    svtype = "DUP"
                    alt = "<DUP>"
                    svlen = end - start + 1
                    if reads_ratio > 1.8:
                        cn = 4
                        gt = "1/1"
                    else:
                        cn = 3
                        gt = "0/1"
                else:
                    logging.warning(f"Unknown CNV type '{call_type}' in row: {row}")
                    continue

                filt = "PASS" if bf >= min_bf else "LowQual"
                qual = f"{min(1000.0, max(0.0, bf)):.2f}"
                var_id = f"EXOMEDEPTH_{sample_id}_{chrom}_{start}_{end}_{svtype}"

                info = (
                    f"SVTYPE={svtype};END={end};SVLEN={svlen};TOOL=EXOMEDEPTH;"
                    f"BF={bf:.2f};READS_RATIO={reads_ratio:.3f};"
                    f"READS_EXPECTED={reads_exp:.1f};READS_OBSERVED={reads_obs:.1f};"
                    f"NEXONS={nexons}"
                )
                fmt = "GT:BF:RR:CN"
                sample_val = f"{gt}:{bf:.2f}:{reads_ratio:.3f}:{cn}"

                records.append((chrom, start, end, var_id, alt, qual, filt, info, fmt, sample_val))

    # Sort records by coordinate: natural chromosome sort, then start position
    def chrom_key(c):
        clean_c = c.lower().replace("chr", "")
        if clean_c.isdigit():
            return (0, int(clean_c))
        elif clean_c == "x":
            return (1, 23)
        elif clean_c == "y":
            return (1, 24)
        elif clean_c in ["m", "mt"]:
            return (1, 25)
        return (2, clean_c)

    records.sort(key=lambda r: (chrom_key(r[0]), r[1]))

    with open(out_vcf_path, 'w') as vcf_out:
        write_vcf_header(vcf_out, contigs, sample_id)
        for chrom, start, end, var_id, alt, qual, filt, info, fmt, sample_val in records:
            vcf_out.write(f"{chrom}\t{start}\t{var_id}\tN\t{alt}\t{qual}\t{filt}\t{info}\t{fmt}\t{sample_val}\n")

    print(f"[INFO] Wrote {len(records)} ExomeDepth variants to {out_vcf_path}")
    return out_vcf_path


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = parse_args()
    convert_csv_to_vcf(
        input_csv=args.input_file,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
        fai_file=args.fai_file,
        min_bf=args.min_bf
    )


if __name__ == "__main__":
    main()
