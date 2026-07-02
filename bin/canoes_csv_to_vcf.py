#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: canoes_csv_to_vcf.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

import argparse  # For parsing command-line arguments
import os  # For file and directory operations
import csv  # For reading and writing CSV files
import logging  # For logging errors and other messages
from datetime import datetime  # For getting the current date and time


def setup_logging(log_file=None):
    """Set up logging configuration."""
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
    """Read sample IDs from a .txt file."""
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
    """Safely convert a value to an integer."""
    try:
        return int(value)
    except ValueError:
        return 0


def safe_float(value):
    """Safely convert a value to a float."""
    try:
        return float(value)
    except ValueError:
        return 0.0


def convert_canoes_csv_to_dict(input_file, sample_file, log_file=None):
    """Convert CANOES CSV to a dictionary of mutations organized by sample ID."""
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
    """Extract reference name from FASTA index file."""
    base_name = os.path.basename(fai_file)
    ref_name = base_name.rsplit('.', 2)[0]
    return ref_name


def create_vcf_contig_lines(fai_file):
    """Create VCF contig lines from FASTA index file."""
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
    """Write the VCF header to the file."""
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


def write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample, log_file=None):
    """Write mutations to the VCF file."""
    if log_file:
        setup_logging(log_file)

    mutation_list, mutations_by_sample = convert_canoes_csv_to_dict(input_file, sample_file, log_file)

    if individual_sample not in mutations_by_sample:
        return

    for mutation in mutations_by_sample[individual_sample]:
        chr_num = mutation['CHR']
        interval = mutation['INTERVAL']

        ref = 'N'
        alt = '<DEL>,<DUP>'

        try:
            start = interval.split(':')[1].split('-')[0]
            end = interval.split(':')[1].split('-')[1]
            svlen = safe_int(end) - safe_int(start) + 1
        except (IndexError, ValueError):
            logging.error(f"Malformed INTERVAL '{interval}' for sample {individual_sample}; skipping record")
            continue

        info = (
            f"END={end};SVLEN={svlen};SVTYPE=CNV;TOOL=CANOES;"
            f"TARGETS={mutation['TARGETS']};NUM_TARG={mutation['NUM_TARG']};MLCN={mutation['MLCN']}"
        )
        filter_status = 'PASS' if safe_float(mutation['Q_SOME']) >= 80 else 'LowQuality'
        
        gt_val = "0"
        if mutation['CNV'] == "DEL":
            gt_val = "1"
        elif mutation['CNV'] == "DUP":
            gt_val = "2"
            
        vcf_file.write(f"chr{chr_num}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\tGT:Q_SOME\t{gt_val}:{safe_float(mutation['Q_SOME'])}\n")


def process_canoes_data(input_file, sample_file, fai_file, output_dir, log_file=None):
    """Process CANOES data and generate individual VCF files for each sample."""
    if log_file:
        setup_logging(log_file)

    sample_list = read_sample_list(sample_file)
    mutation_list, mutations_by_sample = convert_canoes_csv_to_dict(input_file, sample_file, log_file)

    if mutation_list:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for sample in sample_list:
            if not mutations_by_sample[sample]:
                continue
            individual_vcf_path = os.path.join(output_dir, f"{sample}_CANOES_output.vcf")
            with open(individual_vcf_path, 'w') as vcf_file:
                write_vcf_header(vcf_file, fai_file, sample)
                write_vcf_mutations(vcf_file, input_file, sample_file, individual_sample=sample, log_file=None)


def main():
    """Main function to process CANOES data and generate VCF files."""
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
