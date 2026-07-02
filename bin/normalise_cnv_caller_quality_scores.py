#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: normalise_cnv_caller_quality_scores.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

import argparse
import pysam


def standardize_cnv_qual(input_vcf, output_vcf, caller):
    """
    Standardizes CNV quality scores across 7 callers for Truvari integration.
    Valid callers: 'CANOES', 'CLAMMS', 'XHMM', 'GATK', 'DRAGEN', 'CNVKIT', 'INDELIBLE'
    """
    # Map caller names to standard uppercase keys to support aliases and case-insensitivity
    caller = caller.upper()
    if caller in ['GATK-GCNV', 'GCNV']:
        caller = 'GATK'
    elif caller in ['CNV-KIT']:
        caller = 'CNVKIT'
    elif caller in ['DRAGEN_GERMLINE', 'DRAGEN-GERMLINE']:
        caller = 'DRAGEN'

    vcf = pysam.VariantFile(input_vcf, "r")

    # Add new FORMAT headers to preserve original metrics
    vcf.header.formats.add("OQ", "1", "Float", "Original Quality Score")
    vcf.header.formats.add("OAS", "1", "String", "Original Algorithm Score Metric Used")

    out = pysam.VariantFile(output_vcf, "w", header=vcf.header)
    universal_baseline = 100.0

    for record in vcf:
        # Clean up record ID if it contains colons to prevent SURVIVOR merge format field parse error
        if record.id and ":" in record.id:
            record.id = record.id.replace(":", "_")

        # Get the first sample's FORMAT fields for single-sample VCFs
        sample = next(iter(record.samples.values()), None)

        # Safely extract existing QUAL
        orig_qual = record.qual
        if orig_qual is not None and sample is not None:
            sample["OQ"] = orig_qual  # Move to OQ field
            record.qual = None  # Clear the QUAL field to prepare for normalized score

        qual_norm = 0.0
        metric_used = "UNKNOWN"

        try:
            if caller == "CANOES":
                # Q_SOME is a FORMAT field in CANOES VCFs
                if sample is not None and "Q_SOME" in sample:
                    q_some = float(sample["Q_SOME"])
                    qual_norm = min(1000.0, universal_baseline * (q_some / 80.0))
                    metric_used = "Q_SOME"

            elif caller == "CLAMMS":
                # Q_SOME and Q_EXACT are FORMAT fields in CLAMMS VCFs
                if sample is not None and "Q_SOME" in sample and "Q_EXACT" in sample:
                    q_exact = float(sample["Q_EXACT"])
                    q_some = float(sample["Q_SOME"])
                    if q_exact >= 0.0:
                        qual_norm = min(1000.0, universal_baseline * (q_some / 500.0))
                    metric_used = "Q_SOME"

            elif caller == "XHMM":
                if sample is not None and "SQ" in sample and "EQ" in sample and "NDQ" in sample:
                    # GT is a tuple, e.g. (1,) or (2,). Handle index based on genotype:
                    # index 0 for DEL (GT=1), index 1 for DUP (GT=2)
                    gt = sample.get("GT", (None,))
                    idx = 0
                    if gt and any(g == 2 for g in gt if g is not None):
                        idx = 1
                    
                    sq_val = sample["SQ"]
                    eq_val = sample["EQ"]
                    ndq_val = sample["NDQ"]
                    
                    if sq_val is not None and eq_val is not None and ndq_val is not None:
                        sq = float(sq_val[idx]) if isinstance(sq_val, (tuple, list)) else float(sq_val)
                        eq = float(eq_val[idx]) if isinstance(eq_val, (tuple, list)) else float(eq_val)
                        ndq = float(ndq_val[idx]) if isinstance(ndq_val, (tuple, list)) else float(ndq_val)
                        
                        if eq >= 60.0 and ndq >= 60.0:
                            qual_norm = min(1000.0, universal_baseline * (sq / 60.0))
                        metric_used = "SQ"

            elif caller == "GATK":
                if sample is not None and "QS" in sample and "CN" in sample and "NP" in sample:
                    qs = float(sample["QS"])
                    cn = int(sample["CN"])
                    n_int = float(sample["NP"])

                    t_gatk = None
                    if cn == 0:
                        t_gatk = min(1000.0, max(400.0, 10.0 * n_int))
                    elif cn == 1:
                        t_gatk = min(1000.0, max(100.0, 10.0 * n_int))
                    elif cn > 2:
                        t_gatk = min(400.0, max(50.0, 4.0 * n_int))

                    if t_gatk is not None:
                        qual_norm = min(1000.0, universal_baseline * (qs / t_gatk))
                        metric_used = "QS"

            elif caller == "CNVKIT":
                if sample is not None:
                    if "CNQ" in sample:
                        cnq = float(sample["CNQ"])
                        qual_norm = min(1000.0, universal_baseline * (cnq / 5.0))
                        metric_used = "CNQ"
                    elif "GQ" in sample:
                        gq = float(sample["GQ"])
                        qual_norm = min(1000.0, universal_baseline * (gq / 5.0))
                        metric_used = "GQ"

            elif caller == "DRAGEN":
                if orig_qual is not None and "PASS" in record.filter.keys():
                    dragen_qual = float(orig_qual)
                    if dragen_qual <= 10.0:
                        qual_norm = dragen_qual * 10.0
                    else:
                        qual_norm = 100.0 + (dragen_qual - 10.0) * (900.0 / 190.0)
                    metric_used = "QUAL"

            elif caller == "INDELIBLE":
                # SR_TOTAL and AVG_MAPQ are INFO fields in INDELIBLE VCFs (uppercase).
                # MUM_SR and DAD_SR are optional de-novo INFO fields; default to 0
                # when the TSV pre-dates the de-novo analysis step.
                if "SR_TOTAL" in record.info and "AVG_MAPQ" in record.info:
                    sr_total = float(record.info["SR_TOTAL"])
                    avg_mapq = float(record.info["AVG_MAPQ"])
                    mum_sr = float(record.info["MUM_SR"]) if "MUM_SR" in record.info else 0.0
                    dad_sr = float(record.info["DAD_SR"]) if "DAD_SR" in record.info else 0.0

                    if sr_total >= 5.0 and avg_mapq >= 20.0 and mum_sr < 2.0 and dad_sr < 2.0:
                        synthetic_score = sr_total * (avg_mapq / 60.0) * 100.0
                        qual_norm = min(1000.0, universal_baseline * (synthetic_score / 16.67))
                        metric_used = "SYNTHETIC_PHRED"

        except (ValueError, TypeError) as e:
            print(f"Error processing record {record.id}: {e}")  # Logging error

        # Set the normalized quality score in the QUAL field
        record.qual = round(qual_norm, 2)
        if sample is not None:
            sample["OAS"] = metric_used  # Store metric used

        out.write(record)

    vcf.close()
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Normalise CNV caller quality scores to a common scale.'
    )
    parser.add_argument(
        '--input_vcf', required=True,
        help='Input VCF file (plain or bgzipped)'
    )
    parser.add_argument(
        '--output_vcf', required=True,
        help='Output VCF file path'
    )
    parser.add_argument(
        '--caller', required=True,
        help='Name of the CNV caller that produced the input VCF'
    )
    args = parser.parse_args()
    standardize_cnv_qual(args.input_vcf, args.output_vcf, args.caller)
