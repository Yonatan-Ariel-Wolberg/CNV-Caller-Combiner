#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: clamms_bed_to_vcf.py
# Purpose: This pipeline helper script converts CLAMMS caller output BED files
#   into individual-sample VCF files containing called CNV mutations.
# Inputs:
#   - --input_file: The tab-separated CLAMMS output BED file.
#   - --sample_file: A text file listing sample IDs, one per line.
#   - --output_dir: Target directory path to output the generated VCFs.
#   - --fai_file: FASTA index (.fai) reference file to construct VCF headers.
#   - --log_file: Optional path to output error logs.
# Outputs:
#   - Individual VCF files named <sample_id>_CLAMMS_output.vcf in the output directory.
# Operational notes:
#   - Reads and maps BED coordinates to standard VCF format (POS start is extracted).
#   - Maps standard ALT alleles to <DEL> and <DUP> representing deletion and duplication events.
#   - Translates genotype information (1 for DEL, 2 for DUP, 0 fallback) and extracts quality parameters
#     (Q_SOME, Q_EXACT).
#   - Builds chromosome contig headers from the FAI index for canonical chromosomes (chr1-chr22).
#   - Filters variants using the standard CLAMMS quality rules: PASS if Q_SOME >= 500 and Q_EXACT >= 0,
#     otherwise LowQuality.
# =============================================================================

# Import necessary libraries
import argparse  # For parsing command-line arguments
import os  # For file and directory operations
import csv  # For reading and writing CSV files
import logging  # For logging errors and other messages
from datetime import datetime  # For getting the current date and time


def setup_logging(log_file=None):
    """
    Set up standard logging configuration.

    Parameters:
    - log_file: Optional log file path to record errors.
    """
    if log_file:
        logging.basicConfig(
            filename=log_file,
            level=logging.ERROR,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
    else:
        logging.basicConfig(
            level=logging.ERROR,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
    logging.getLogger().addHandler(logging.StreamHandler())


def read_sample_list(sample_file):
    """
    Read unique sample IDs from a text file, one sample per line.

    Parameters:
    - sample_file: Path to the text file containing sample IDs.

    Returns:
    - List of unique, stripped sample ID strings.
    """
    sample_list = []
    try:
        with open(sample_file, 'r') as f:
            for line in f:
                # Read sample IDs, strip whitespace and maintain order
                sample_id = line.strip()
                if sample_id and sample_id not in sample_list:
                    sample_list.append(sample_id)

    except IOError as e:
        logging.error(f"Error reading sample IDs from {sample_file}: {e}")

    return sample_list


def safe_int(value):
    """
    Safely convert a string/float value to an integer, returning 0 on failure.

    Parameters:
    - value: The input value to convert.

    Returns:
    - Converted integer value or 0.
    """
    try:
        return int(value)
    except ValueError:
        return 0


def safe_float(value):
    """
    Safely convert a string value to a float, returning 0.0 on failure.

    Parameters:
    - value: The input value to convert.

    Returns:
    - Converted float value or 0.0.
    """
    try:
        return float(value)
    except ValueError:
        return 0.0


def convert_clamms_bed_to_dict(input_file, sample_file, log_file=None):
    """
    Convert a CLAMMS BED file to a dictionary of mutations organized by sample ID.

    The BED file is expected to be tab-separated, containing fields:
      0: CHROM, 1: START, 2: END, 3: INTERVAL, 4: SAMPLE, 5: CNV_TYPE, 6: MLCN,
      7: NUM_WINDOWS, 8: Q_SOME, 9: Q_EXACT, 10: Q_LEFT_EXTEND, 11: LEFT_EXTEND_COORD,
      12: Q_RIGHT_EXTEND, 13: RIGHT_EXTEND_COORD, 14: Q_LEFT_CONTRACT, 15: LEFT_CONTRACT_COORD,
      16: Q_RIGHT_CONTRACT, 17: RIGHT_CONTRACT_COORD.

    Parameters:
    - input_file: Path to the input BED file.
    - sample_file: Path to the sample list file.
    - log_file: Optional log file to record errors.

    Returns:
    - A tuple containing:
      1. List of all mutations parsed.
      2. Dict mapping sample IDs to lists of their respective mutations.
    """
    # Set up logging
    if log_file:
        setup_logging(log_file)

    sample_list = read_sample_list(sample_file)

    mutations_by_sample = {sample: [] for sample in sample_list}
    mutation_list = []

    try:
        with open(input_file, 'r') as bedfile:
            for row_number, row in enumerate(bedfile):
                fields = row.strip().split('\t')
                if len(fields) < 18:
                    logging.error(f"Row {row_number+1} does not have enough fields: {row.strip()}")
                    continue

                mutation = {
                    'CHROM': fields[0],
                    'START': fields[1],
                    'END': fields[2],
                    'INTERVAL': fields[3],
                    'CNV_TYPE': fields[5],
                    'MLCN': fields[6],
                    'NUM_WINDOWS': fields[7],
                    'Q_SOME': fields[8],
                    'Q_EXACT': fields[9],
                    'Q_LEFT_EXTEND': fields[10],
                    'LEFT_EXTEND_COORD': fields[11],
                    'Q_RIGHT_EXTEND': fields[12],
                    'RIGHT_EXTEND_COORD': fields[13],
                    'Q_LEFT_CONTRACT': fields[14],
                    'LEFT_CONTRACT_COORD': fields[15],
                    'Q_RIGHT_CONTRACT': fields[16],
                    'RIGHT_CONTRACT_COORD': fields[17],
                }

                mutation_list.append(mutation)

                sample_id = fields[4]
                if sample_id in mutations_by_sample:
                    mutations_by_sample[sample_id].append(mutation)
                else:
                    logging.warning(
                        f"Sample '{sample_id}' found in BED but missing from sample list; "
                        "mutation will be ignored"
                    )

    except IOError as e:
        logging.error(f"Error reading input BED file {input_file}: {e}")
        return [], {}

    return mutation_list, mutations_by_sample


def extract_ref_name(fai_file):
    """
    Extract reference name from the filename of the FASTA index file.

    Parameters:
    - fai_file: Path to the FASTA index file.

    Returns:
    - String representing reference genome name.
    """
    base_name = os.path.basename(fai_file)
    ref_name = base_name.rsplit('.', 2)[0]
    return ref_name


def create_vcf_contig_lines(fai_file):
    """
    Create VCF contig definition lines parsed from the FASTA index (.fai) file.
    Filters to keep only canonical chromosomes (chr1-chr22).

    Parameters:
    - fai_file: Path to FASTA index file.

    Returns:
    - List of VCF contig line headers.
    """
    contig_lines = []
    valid_contigs = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'
    ]

    valid_contig_set = set(valid_contigs)

    with open(fai_file) as fai:
        for line in fai:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                logging.error(f"Malformed FAI line (expected >=2 columns): {line.strip()}")
                continue
            contig_name = parts[0]
            try:
                contig_length = int(parts[1])
            except ValueError:
                logging.error(f"Invalid contig length in FAI line: {line.strip()}")
                continue

            if contig_name in valid_contig_set:
                contig_lines.append(f"##contig=<ID={contig_name},length={contig_length}>")

    return contig_lines


def write_vcf_header(vcf_file, fai_file, individual_sample):
    """
    Write complete VCF headers including format attributes, custom info tags, and chromosomes.

    Parameters:
    - vcf_file: VCF file context to write the header into.
    - fai_file: Reference genome FASTA index file.
    - individual_sample: Sample name for the individual VCF.
    """
    ref_name = extract_ref_name(fai_file)
    sorted_contig_lines = create_vcf_contig_lines(fai_file)

    vcf_file.write("##fileformat=VCFv4.1\n")
    vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
    vcf_file.write("##source=CLAMMS\n")
    vcf_file.write(f"##reference={ref_name}\n")
    vcf_file.write("##phasing=partial\n")

    for line in sorted_contig_lines:
        vcf_file.write(line + "\n")

    vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    vcf_file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    vcf_file.write('##FILTER=<ID=LowQuality,Description="Low quality SV">\n')
    vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
    vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    vcf_file.write("##INFO=<ID=TOOL,Number=1,Type=String,Description=\"CNV Caller used to call CNV\">\n")
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
    vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strand direction\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Quality score of the CNV event\">\n")
    vcf_file.write("##FORMAT=<ID=Q_EXACT,Number=1,Type=Float,Description=\"Exact boundary quality score of the CNV event\">\n")

    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + individual_sample + "\n")


def write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample, log_file=None, sample_mutations=None):
    """
    Write variant records for an individual sample to the VCF file.

    Parameters:
    - vcf_file: File handle of the output VCF.
    - input_file: Path to CLAMMS BED output.
    - sample_file: Path to sample list file.
    - individual_sample: Sample ID string.
    - log_file: Optional log path.
    - sample_mutations: Optional pre-parsed mutations list for individual_sample.
    """
    if log_file:
        setup_logging(log_file)

    if sample_mutations is None:
        mutation_list, mutations_by_sample = convert_clamms_bed_to_dict(input_file, sample_file, log_file)
        sample_mutations = mutations_by_sample.get(individual_sample, [])

    for mutation in sample_mutations:
        chr_num = mutation['CHROM']

        cnv = mutation['CNV_TYPE']
        start = mutation['START']
        end = mutation['END']
        mlcn = mutation['MLCN']
        num_windows = mutation['NUM_WINDOWS']
        q_some = safe_float(mutation['Q_SOME'])
        q_exact = safe_float(mutation['Q_EXACT'])

        svmethod = 'CLAMMS'
        strand = '.'
        ref = 'N'
        alt = f"<{cnv}>" if cnv in ["DEL", "DUP"] else '<DEL>,<DUP>'

        svlen = safe_int(end) - safe_int(start) + 1

        info = f"END={end};SVLEN={svlen};SVTYPE={cnv};TOOL={svmethod};STRANDS={strand}"

        # Standard CLAMMS filters: PASS if Q_SOME >= 50 and Q_EXACT >= 0 (or >= -0.5 for DUP), otherwise LowQuality
        filter_status = 'PASS' if q_some >= 50 and ((cnv == 'DEL' and q_exact >= 0) or (cnv == 'DUP' and q_exact >= -0.5) or (cnv not in ['DEL', 'DUP'] and q_exact >= 0)) else 'LowQuality'
        format_field = "GT:Q_SOME:Q_EXACT"
        chrom_formatted = f"chr{chr_num}" if not str(chr_num).startswith("chr") else chr_num

        gt_val = "1" if cnv == "DEL" else ("2" if cnv == "DUP" else "0")

        vcf_file.write(f"{chrom_formatted}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\t{format_field}\t{gt_val}:{q_some}:{q_exact}\n")


def process_clamms_data(input_file, sample_file, fai_file, output_dir, log_file=None):
    """
    Load data for all samples, and generate individual output VCF files for each.

    Parameters:
    - input_file: Tabular input file path.
    - sample_file: Target samples list.
    - fai_file: Reference FAI index.
    - output_dir: Target output folder path.
    - log_file: Optional logger path.
    """
    if log_file:
        setup_logging(log_file)

    sample_list = read_sample_list(sample_file)
    mutation_list, mutations_by_sample = convert_clamms_bed_to_dict(input_file, sample_file, log_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for sample in sample_list:
        sample_mutations = mutations_by_sample.get(sample, [])
        if not sample_mutations:
            continue
        individual_vcf_path = os.path.join(output_dir, f"{sample}_CLAMMS_output.vcf")
        with open(individual_vcf_path, 'w') as vcf_file:
            write_vcf_header(vcf_file, fai_file, individual_sample=sample)
            write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample=sample, log_file=None, sample_mutations=sample_mutations)


def main():
    """Main execution entry point parser."""
    parser = argparse.ArgumentParser(description="Process CLAMMS BED to VCF conversion.")
    parser.add_argument('--input_file', type=str, required=True, help="Path to CLAMMS BED file.")
    parser.add_argument('--sample_file', type=str, required=True, help="Path to sample list file.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save output VCF files.")
    parser.add_argument('--fai_file', type=str, required=True, help="FASTA Index File")
    parser.add_argument('--log_file', type=str, default=None, help="Optional log file.")

    args = parser.parse_args()

    process_clamms_data(args.input_file, args.sample_file, args.fai_file, args.output_dir, log_file=args.log_file)


if __name__ == "__main__":
    main()
