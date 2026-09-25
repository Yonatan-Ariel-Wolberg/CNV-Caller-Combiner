#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: canoes_csv_to_vcf.py
# Purpose: This pipeline helper script converts CANOES caller tabular CSV outputs
#   into standard individual-sample VCF files containing called CNV mutations.
# Inputs:
#   - --input_file: Tab-separated CANOES CSV output containing called variants.
#   - --sample_file: A text file containing a list of sample IDs, one per line.
#   - --output_dir: Target directory path to write the generated VCF files.
#   - --fai_file: FASTA index (.fai) reference file to construct VCF contig headers.
#   - --log_file: Optional path to write runtime error logs.
# Outputs:
#   - Individual VCF files named <sample_id>_CANOES_output.vcf in the output directory.
# Operational notes:
#   - Parses genomic coordinate intervals from the "INTERVAL" CSV field.
#   - Maps standard ALT alleles to <DEL> and <DUP> representing deletion and duplication events.
#   - Infers Genotype (GT) field value (1 for DEL, 2 for DUP, 0 fallback) and extracts quality scores (Q_SOME).
#   - Automatically extracts reference names and builds the required contig headers from the FAI index.
#   - Applies a default Q_SOME filter (PASS if Q_SOME >= 80, otherwise LowQuality).
# =============================================================================

import argparse  # For parsing command-line arguments
import os  # For file and directory operations
import csv  # For reading and writing CSV files
import logging  # For logging errors and other messages
from datetime import datetime  # For getting the current date and time


def setup_logging(log_file=None):
    """
    Set up standard logging configuration.

    Parameters:
    - log_file: Path to write logging output. If None, prints to stream.
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
                sample_id = line.strip()
                if sample_id and sample_id not in sample_list:
                    sample_list.append(sample_id)
    except IOError as e:
        logging.error(f"Error reading sample IDs from {sample_file}: {e}")
        return []

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


def convert_canoes_csv_to_dict(input_file, sample_file, log_file=None):
    """
    Parse a CANOES CSV file and categorize mutation records by their Sample ID.

    The CSV contains tab-separated columns:
      SAMPLE, CNV, INTERVAL, KB, CHR, MID_BP, TARGETS, NUM_TARG, MLCN, Q_SOME.
      
    Args:
        input_file (str): Path to the CANOES output CSV file.
        sample_file (str): Path to the sample list file.
        log_file (str, optional): Optional path to output error logs.

    Returns:
    - A tuple containing:
      1. List of all unique mutations processed.
      2. Dict mapping sample IDs to lists of their respective mutations.
    """
    # Set up logging
    if log_file:
        setup_logging(log_file)

    sample_list = read_sample_list(sample_file)
    mutations_by_sample = {sample: [] for sample in sample_list}
    mutation_list = []

    try:
        with open(input_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                if row['SAMPLE'] in sample_list:
                    mutation = {
                        'CNV': row['CNV'],
                        'INTERVAL': row['INTERVAL'],
                        'KB': row['KB'],
                        'CHR': row['CHR'],
                        'MID_BP': row['MID_BP'],
                        'TARGETS': row['TARGETS'],
                        'NUM_TARG': row['NUM_TARG'],
                        'SAMPLE': row['SAMPLE'],
                        'MLCN': row['MLCN'],
                        'Q_SOME': row['Q_SOME']
                    }

                    if mutation not in mutation_list:
                        mutation_list.append(mutation)

                    mutations_by_sample[row['SAMPLE']].append(mutation)

    except IOError as e:
        logging.error(f"Error reading input CSV file {input_file}: {e}")
        return [], {}
    except Exception as e:
        logging.error(f"Error processing input CSV file {input_file}: {e}")
        return [], {}

    return mutation_list, mutations_by_sample


def extract_ref_name(fai_file):
    """
    Extract a clean reference name string based on the filename of the FAI file.
    
    Example:
      /path/to/hg19.fa.fai -> hg19

    Parameters:
    - fai_file: Path to the FASTA index file (e.g. hg19.fa.fai).

    Returns:
    - String representing reference genome name (e.g. hg19).
    """
    base_name = os.path.basename(fai_file)
    ref_name = base_name.rsplit('.', 2)[0]
    return ref_name


def create_vcf_contig_lines(fai_file):
    """
    Create VCF contig definition lines parsed from the FASTA index (.fai) file.
    
    Reads lines of FAI, which contain tab-separated columns:
      <contig_name>\t<length>\t...
    Generates VCF contig lines in the format: ##contig=<ID=contig_name,length=contig_length>

    Parameters:
    - fai_file: Path to the FASTA index file.

    Returns:
    - List of VCF header contig strings.
    """
    contig_lines = []

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

            contig_line = f"##contig=<ID={contig_name},length={contig_length}>"
            contig_lines.append(contig_line)

    return contig_lines


def write_vcf_header(vcf_file, fai_file, sample_id):
    """
    Write VCF headers including formats, custom tags, reference, and contigs.

    Parameters:
    - vcf_file: File handle of the output VCF.
    - fai_file: Path to the reference index file.
    - sample_id: Name of the individual sample.
    """
    ref_name = extract_ref_name(fai_file)
    sorted_contig_lines = create_vcf_contig_lines(fai_file)

    vcf_file.write("##fileformat=VCFv4.1\n")
    vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
    vcf_file.write("##source=CANOES\n")
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
    vcf_file.write("##INFO=<ID=TARGETS,Number=1,Type=String,Description=\"Probes targeted by the SV\">\n")
    vcf_file.write("##INFO=<ID=NUM_TARG,Number=1,Type=Integer,Description=\"Number of targets involved\">\n")
    vcf_file.write("##INFO=<ID=MLCN,Number=1,Type=Integer,Description=\"Most likely copy number\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Quality score of the CNV event\">\n")

    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")


def write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample, log_file=None, sample_mutations=None):
    """
    Translate individual mutations into standard VCF records and write to file.
    
    - Parses INTERVAL (e.g. '1:10000-20000') into POS and END.
    - Calculates SVLEN = END - START + 1.
    - Filters: PASS if Q_SOME >= 80, otherwise LowQuality.
    - Genotypes: 
        * DEL maps to GT=1
        * DUP maps to GT=2
        * Other maps to GT=0

    Parameters:
    - vcf_file: Output file handle to write records into.
    - input_file: Path to CANOES CSV input.
    - sample_file: Path to sample list file.
    - individual_sample: Target sample ID being written.
    - log_file: Optional logging output path.
    - sample_mutations: Optional pre-parsed mutations list for individual_sample.
    """
    if log_file:
        setup_logging(log_file)

    if sample_mutations is None:
        mutation_list, mutations_by_sample = convert_canoes_csv_to_dict(input_file, sample_file, log_file)
        sample_mutations = mutations_by_sample.get(individual_sample, [])

    for mutation in sample_mutations:
        chr_num = str(mutation['CHR'])
        chrom_formatted = chr_num if chr_num.startswith("chr") else f"chr{chr_num}"
        interval = mutation['INTERVAL']

        ref = 'N'
        alt = f"<{mutation['CNV']}>" if mutation['CNV'] in ["DEL", "DUP"] else '<DEL>,<DUP>'

        # Split coordinate details from INTERVAL column (format expected: chromosome:start-end)
        try:
            start = interval.split(':')[1].split('-')[0]
            end = interval.split(':')[1].split('-')[1]
            svlen = safe_int(end) - safe_int(start) + 1
        except (IndexError, ValueError):
            logging.error(f"Malformed INTERVAL '{interval}' for sample {individual_sample}; skipping record")
            continue

        # Format custom INFO, filter flags, and format parameters
        info = (
            f"END={end};SVLEN={svlen};SVTYPE=CNV;TOOL=CANOES;"
            f"TARGETS={mutation['TARGETS']};NUM_TARG={mutation['NUM_TARG']};MLCN={mutation['MLCN']}"
        )
        
        # Apply standard quality filtration threshold (Q_SOME >= 80.0 passes)
        filter_status = 'PASS' if safe_float(mutation['Q_SOME']) >= 80 else 'LowQuality'
        
        # Determine genotype parameter
        # GT represents index of ALT list: since it is biallelic, GT=1 points to the CNV call
        gt_val = "1" if mutation['CNV'] in ["DEL", "DUP"] else "0"
            
        vcf_file.write(f"{chrom_formatted}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\tGT:Q_SOME\t{gt_val}:{safe_float(mutation['Q_SOME'])}\n")


def process_canoes_data(input_file, sample_file, fai_file, output_dir, log_file=None):
    """
    Iterate over the sample list and coordinate data to write separate VCF outputs.

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
    mutation_list, mutations_by_sample = convert_canoes_csv_to_dict(input_file, sample_file, log_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Generate separate VCFs for each sample in the project list
    for sample in sample_list:
        if not mutations_by_sample.get(sample):
            continue
        individual_vcf_path = os.path.join(output_dir, f"{sample}_CANOES_output.vcf")
        with open(individual_vcf_path, 'w') as vcf_file:
            write_vcf_header(vcf_file, fai_file, sample)
            write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample=sample, log_file=None, sample_mutations=mutations_by_sample.get(sample, []))


def main():
    """Main execution block parsing command-line parameters."""
    parser = argparse.ArgumentParser(description="Process CANOES CSV to VCF conversion.")
    parser.add_argument('--input_file', type=str, help="Path to CANOES CSV file.")
    parser.add_argument('--sample_file', type=str, help="Path to sample list file.")
    parser.add_argument('--output_dir', type=str, help="Directory to save output VCF files.")
    parser.add_argument('--fai_file', type=str, required=True, help="FASTA Index File")
    parser.add_argument('--log_file', type=str, default=None, help="Optional log file.")

    args = parser.parse_args()

    process_canoes_data(args.input_file, args.sample_file, args.fai_file, args.output_dir, log_file=args.log_file)


if __name__ == "__main__":
    main()
