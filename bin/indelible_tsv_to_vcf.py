#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: indelible_tsv_to_vcf.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

import argparse
import os
import csv
import logging
from datetime import datetime


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


def extract_ref_name(fai_file):
    """Extract reference name from FASTA index file."""
    if not fai_file or not os.path.exists(fai_file):
        return "unknown_reference"
    
    base_name = os.path.basename(fai_file)
    ref_name = base_name.rsplit('.', 2)[0]
    return ref_name


def create_vcf_contig_lines(fai_file):
    """Create VCF contig lines from FASTA index file."""
    if not fai_file or not os.path.exists(fai_file):
        return []

    contig_lines = []
    valid_contigs = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
        'chrX', 'chrY'
    ]
    valid_contig_set = set(valid_contigs)

    with open(fai_file) as fai:
        for line in fai:
            parts = line.strip().split("\t")
            contig_name = parts[0]
            contig_length = int(parts[1])

            if contig_name in valid_contig_set:
                contig_line = f"##contig=<ID={contig_name},length={contig_length}>"
                contig_lines.append((valid_contigs.index(contig_name), contig_line))

    contig_lines.sort()
    return [line for _, line in contig_lines]


def write_vcf_header(vcf_file, fai_file, sample_id):
    """Write the VCF header to the file."""
    ref_name = extract_ref_name(fai_file)
    sorted_contig_lines = create_vcf_contig_lines(fai_file)

    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
    vcf_file.write("##source=INDELIBLE\n")
    vcf_file.write(f"##reference={ref_name}\n")
    vcf_file.write("##phasing=partial\n")

    for line in sorted_contig_lines:
        vcf_file.write(line + "\n")

    vcf_file.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
    vcf_file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    vcf_file.write('##FILTER=<ID=LowQuality,Description="Low quality SV">\n')
    vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
    vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    vcf_file.write("##INFO=<ID=TOOL,Number=1,Type=String,Description=\"CNV Caller used to call CNV\">\n")
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
    vcf_file.write("##INFO=<ID=SR_TOTAL,Number=1,Type=Integer,Description=\"Total split reads supporting the variant\">\n")
    vcf_file.write("##INFO=<ID=COVERAGE,Number=1,Type=Integer,Description=\"Coverage at the position\">\n")
    vcf_file.write("##INFO=<ID=PROB_Y,Number=1,Type=Float,Description=\"Probability score from INDELIBLE\">\n")
    vcf_file.write("##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Inserted Sequence\">\n")
    vcf_file.write("##INFO=<ID=INSERTION_CONTEXT,Number=1,Type=String,Description=\"Insertion context\">\n")
    vcf_file.write("##INFO=<ID=DELETION_CONTEXT,Number=1,Type=String,Description=\"Deletion context\">\n")
    vcf_file.write("##INFO=<ID=SR_TOTAL_LONG,Number=1,Type=Integer,Description=\"Total long split reads\">\n")
    vcf_file.write("##INFO=<ID=SR_TOTAL_SHORT,Number=1,Type=Integer,Description=\"Total short split reads\">\n")
    vcf_file.write("##INFO=<ID=SR_LONG_5,Number=1,Type=Integer,Description=\"Long split reads 5 prime\">\n")
    vcf_file.write("##INFO=<ID=SR_SHORT_5,Number=1,Type=Integer,Description=\"Short split reads 5 prime\">\n")
    vcf_file.write("##INFO=<ID=SR_LONG_3,Number=1,Type=Integer,Description=\"Long split reads 3 prime\">\n")
    vcf_file.write("##INFO=<ID=SR_SHORT_3,Number=1,Type=Integer,Description=\"Short split reads 3 prime\">\n")
    vcf_file.write("##INFO=<ID=SR_ENTROPY,Number=1,Type=Float,Description=\"Split read entropy\">\n")
    vcf_file.write("##INFO=<ID=CONTEXT_ENTROPY,Number=1,Type=Float,Description=\"Context entropy\">\n")
    vcf_file.write("##INFO=<ID=ENTROPY_UPSTREAM,Number=1,Type=Float,Description=\"Upstream entropy\">\n")
    vcf_file.write("##INFO=<ID=ENTROPY_DOWNSTREAM,Number=1,Type=Float,Description=\"Downstream entropy\">\n")
    vcf_file.write("##INFO=<ID=SR_SW_SIMILARITY,Number=1,Type=Float,Description=\"Split read Smith-Waterman similarity\">\n")
    vcf_file.write("##INFO=<ID=AVG_AVG_SR_QUAL,Number=1,Type=Float,Description=\"Average split read quality\">\n")
    vcf_file.write("##INFO=<ID=AVG_MAPQ,Number=1,Type=Float,Description=\"Average MAPQ\">\n")
    vcf_file.write("##INFO=<ID=PROB_N,Number=1,Type=Float,Description=\"Probability score N\">\n")
    vcf_file.write("##INFO=<ID=DDG2P,Number=1,Type=String,Description=\"DDG2P gene\">\n")
    vcf_file.write("##INFO=<ID=HGNC,Number=1,Type=String,Description=\"HGNC gene\">\n")
    vcf_file.write("##INFO=<ID=HGNC_CONSTRAINED,Number=1,Type=String,Description=\"HGNC constrained\">\n")
    vcf_file.write("##INFO=<ID=EXONIC,Number=1,Type=String,Description=\"Exonic variant\">\n")
    vcf_file.write("##INFO=<ID=TRANSCRIPTS,Number=1,Type=String,Description=\"Transcripts\">\n")
    vcf_file.write("##INFO=<ID=EXON_NUMBERS,Number=1,Type=String,Description=\"Exon numbers\">\n")
    vcf_file.write("##INFO=<ID=MAF,Number=1,Type=String,Description=\"Minor Allele Frequency\">\n")
    vcf_file.write("##INFO=<ID=BLAST_HIT,Number=1,Type=String,Description=\"Blast hit\">\n")
    vcf_file.write("##INFO=<ID=BLAST_STRAND,Number=1,Type=String,Description=\"Blast strand\">\n")
    vcf_file.write("##INFO=<ID=BLAST_IDENTITY,Number=1,Type=String,Description=\"Blast identity\">\n")
    vcf_file.write("##INFO=<ID=BLAST_DIST,Number=1,Type=String,Description=\"Blast distance\">\n")
    vcf_file.write("##INFO=<ID=BLAST_HGNC,Number=1,Type=String,Description=\"Blast HGNC\">\n")
    vcf_file.write("##INFO=<ID=BLAST_HGNC_CONSTRAINED,Number=1,Type=String,Description=\"Blast HGNC constrained\">\n")
    vcf_file.write("##INFO=<ID=BLAST_DDG2P,Number=1,Type=String,Description=\"Blast DDG2P\">\n")
    vcf_file.write("##INFO=<ID=MUM_SR,Number=1,Type=Integer,Description=\"Split reads in mother supporting the variant (de novo analysis)\">\n")
    vcf_file.write("##INFO=<ID=DAD_SR,Number=1,Type=Integer,Description=\"Split reads in father supporting the variant (de novo analysis)\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_file.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality (Mapped from PROB_Y)\">\n")
    
    vcf_file.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")


def convert_indelible_tsv_to_vcf(input_file, output_dir, sample_id, fai_file, log_file=None):
    """Process INDELIBLE TSV data and generate a VCF file."""
    if log_file:
        setup_logging(log_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_vcf_path = os.path.join(output_dir, f"{sample_id}_INDELIBLE_output.vcf")

    with open(output_vcf_path, 'w') as vcf_file:
        write_vcf_header(vcf_file, fai_file, sample_id)

        try:
            with open(input_file, 'r') as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter='\t')
                
                for row_num, row in enumerate(reader):
                    try:
                        chrom_num = row['chrom']
                        # Strip "chr" if present for uniformity or ensure it's there based on your FAI
                        chrom_formatted = f"chr{chrom_num}" if not str(chrom_num).startswith("chr") else chrom_num
                        
                        pos = row['position']
                        coverage = row.get('coverage', '0')
                        sr_total = row.get('sr_total', '0')
                        seq = row.get('seq_longest', '')
                        predicted = row.get('predicted', 'N')
                        prob_y = safe_float(row.get('prob_Y', '0.0'))

                        # INDELIBLE primary detects insertions
                        svtype = 'INS'
                        ref = 'N'
                        alt = '<INS>'
                        svlen = len(seq) if seq else 0
                        end = pos  # For insertions, END is conventionally the same as POS
                        tool = 'INDELIBLE'

                        # Build the base INFO field
                        info = f"END={end};SVLEN={svlen};SVTYPE={svtype};TOOL={tool};SR_TOTAL={sr_total};COVERAGE={coverage};PROB_Y={prob_y}"
                        if seq and seq != 'NA':
                            info += f";SEQ={seq}"

                        # Dictionary mapping exact VCF INFO tags to TSV column headers
                        optional_fields = {
                            'INSERTION_CONTEXT': 'insertion_context',
                            'DELETION_CONTEXT': 'deletion_context',
                            'SR_TOTAL_LONG': 'sr_total_long',
                            'SR_TOTAL_SHORT': 'sr_total_short',
                            'SR_LONG_5': 'sr_long_5',
                            'SR_SHORT_5': 'sr_short_5',
                            'SR_LONG_3': 'sr_long_3',
                            'SR_SHORT_3': 'sr_short_3',
                            'SR_ENTROPY': 'sr_entropy',
                            'CONTEXT_ENTROPY': 'context_entropy',
                            'ENTROPY_UPSTREAM': 'entropy_upstream',
                            'ENTROPY_DOWNSTREAM': 'entropy_downstream',
                            'SR_SW_SIMILARITY': 'sr_sw_similarity',
                            'AVG_AVG_SR_QUAL': 'avg_avg_sr_qual',
                            'AVG_MAPQ': 'avg_mapq',
                            'PROB_N': 'prob_N',
                            'DDG2P': 'ddg2p',
                            'HGNC': 'hgnc',
                            'HGNC_CONSTRAINED': 'hgnc_constrained',
                            'EXONIC': 'exonic',
                            'TRANSCRIPTS': 'transcripts',
                            'EXON_NUMBERS': 'exon_numbers',
                            'MAF': 'maf',
                            'BLAST_HIT': 'blast_hit',
                            'BLAST_STRAND': 'blast_strand',
                            'BLAST_IDENTITY': 'blast_identity',
                            'BLAST_DIST': 'blast_dist',
                            'BLAST_HGNC': 'blast_hgnc',
                            'BLAST_HGNC_CONSTRAINED': 'blast_hgnc_constrained',
                            'BLAST_DDG2P': 'blast_ddg2p',
                            'MUM_SR': 'mum_sr',
                            'DAD_SR': 'dad_sr'
                        }

                        # Strictly append each optional field to the INFO string (NOT the FORMAT string)
                        for info_key, row_key in optional_fields.items():
                            val = row.get(row_key, '').strip()
                            if val and val != 'NA':
                                # Replace spaces and semicolons to keep VCF parsing robust
                                clean_val = val.replace(' ', '_').replace(';', ',')
                                info += f";{info_key}={clean_val}"

                        # Determine FILTER and Genotype
                        if predicted == 'Y':
                            filter_status = 'PASS'
                            gt = '0/1'
                        else:
                            filter_status = 'LowQuality'
                            gt = '0/0'

                        # The FORMAT field is kept strictly to GT:GQ
                        format_field = "GT:GQ"
                        
                        # Write mutation to the VCF file with fields explicitly separated
                        vcf_file.write(f"{chrom_formatted}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\t{format_field}\t{gt}:{prob_y}\n")
                        
                    except KeyError as e:
                        logging.error(f"Missing expected column at row {row_num}: {e}")
                        
        except Exception as e:
            logging.error(f"Failed to process {input_file}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Convert INDELIBLE TSV to VCF.")
    parser.add_argument('--input_file', type=str, required=True, help="Path to INDELIBLE TSV output file.")
    parser.add_argument('--sample_id', type=str, required=True, help="Sample ID associated with the input file.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save output VCF file.")
    parser.add_argument('--fai_file', type=str, default=None, help="Optional FASTA Index File")
    parser.add_argument('--log_file', type=str, default=None, help="Optional log file.")

    args = parser.parse_args()

    convert_indelible_tsv_to_vcf(
        input_file=args.input_file,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        fai_file=args.fai_file,
        log_file=args.log_file
    )


if __name__ == "__main__":
    main()
