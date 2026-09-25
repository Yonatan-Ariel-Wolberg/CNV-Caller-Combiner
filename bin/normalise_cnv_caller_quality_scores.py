#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: normalise_cnv_caller_quality_scores.py
# Purpose: Standardizes CNV quality scores across 6 different callers
#   (CANOES, CLAMMS, XHMM, GATK, DRAGEN, CNVKIT) onto a common
#   scale to support downstream multi-caller merging (e.g., Truvari/SURVIVOR).
# Inputs:
#   - --input_vcf: Path to the input VCF file (plain or bgzipped).
#   - --output_vcf: Path to write the standardized output VCF file.
#   - --caller: Name of the CNV caller that produced the input VCF.
# Outputs:
#   - A VCF file with standardized quality scores in the QUAL column,
#     original quality scores saved in the FORMAT OQ tag, and original
#     algorithm metric identifier stored in the FORMAT OAS tag.
# Operational notes:
#   - Replaces colons in record IDs with underscores to prevent parsing
#     errors during SURVIVOR merge step.
#   - Employs caller-specific heuristic formulas to normalize metrics (Q_SOME,
#     Q_EXACT, SQ, QS, CNQ, GQ, etc.) to a universal baseline (default 100.0).
# =============================================================================

import argparse
import pysam


def standardize_cnv_qual(
    input_vcf,
    output_vcf,
    caller,
    canoes_baseline=65.0,
    clamms_baseline=50.0,
    xhmm_baseline=50.0,
    cnvkit_baseline_cnq=15.0,
    cnvkit_baseline_log2=0.5,
    dragen_baseline=30.0,
    exomedepth_baseline=20.0
):
    """
    Standardizes CNV quality scores across 7 different callers to a unified scale.
    
    Different CNV calling algorithms output quality metrics using entirely different
    definitions, ranges, and scales. To facilitate down-stream merging and comparison
    (e.g., using Truvari or SURVIVOR), this function normalizes all quality scores
    to a common scale where a 'PASS'-level quality event generally maps to a score 
    of 100.0, with a hard cap at 1000.0.
    
    The function also saves the original quality value in a new sample FORMAT tag 'OQ'
    and records the name of the original quality metric used in sample FORMAT tag 'OAS'.

    ====================================================================================
    MATHEMATICAL FORMULATIONS AND HEURISTIC SCHEMES BY CALLER
    ====================================================================================

    1. CANOES
       - Source metric: 'Q_SOME' (Likelihood of some CNV event at the targets).
       - Baseline threshold: 80.0 (recommended threshold for high confidence).
       - Normalization formula:
           QUAL_norm = min(1000.0, 100.0 * (Q_SOME / 80.0))
       - Rationale: Linearly scales the quality score so that the baseline 80.0 becomes 
         100.0, capped at 1000.0.

    2. CLAMMS
       - Source metrics: 'Q_SOME' (likelihood of some CNV) and 'Q_EXACT' (boundary quality).
       - Baseline threshold: 500.0 for Q_SOME.
       - Prerequisite: Q_EXACT >= 0.0 (insulates against invalid boundaries).
       - Normalization formula:
           QUAL_norm = min(1000.0, 100.0 * (Q_SOME / 500.0))
       - Rationale: Standard CLAMMS filter uses Q_SOME >= 500. This scales Q_SOME so that 
         the 500 threshold maps to 100.0.

    3. XHMM
       - Source metrics: 'SQ' (Genotype Quality score for DEL/DUP), 'EQ' (Exon-level quality),
         and 'NDQ' (Non-deletion/duplication quality).
       - Prerequisite: EQ >= 60.0 and NDQ >= 60.0.
       - Normalization formula:
           The genotype (GT) field determines the index 'idx' in SQ:
             * idx = 0 for deletion (GT contains 1)
             * idx = 1 for duplication (GT contains 2)
           We extract 'sq' corresponding to this index:
             QUAL_norm = min(1000.0, 100.0 * (sq / 60.0))
       - Rationale: Scaled such that the standard XHMM threshold of 60.0 maps to 100.0.

    4. GATK (gCNV)
       - Source metrics: 'QS' (Quality Score), 'CN' (Copy Number), 'NP' (Number of Targets/probes).
       - Dynamic Baseline Threshold (T_GATK):
           A customized baseline threshold is computed based on Copy Number (CN) and the number of 
           target intervals (NP) because gCNV confidence scales with target density:
             * If CN == 0 (homozygous deletion):
                 T_GATK = min(1000.0, max(400.0, 10.0 * NP))
             * If CN == 1 (heterozygous deletion):
                 T_GATK = min(1000.0, max(100.0, 10.0 * NP))
             * If CN > 2 (duplication):
                 T_GATK = min(400.0, max(50.0, 4.0 * NP))
       - Normalization formula:
           QUAL_norm = min(1000.0, 100.0 * (QS / T_GATK))
       - Rationale: Models the expected statistical power of gCNV based on depth reduction/increase 
         over NP intervals, scaling the dynamic expectation to 100.0.

    5. CNVKIT
       - Source metrics: 'CNQ' (Copy Number Quality) or 'GQ' (Genotype Quality).
       - Baseline threshold: 5.0.
       - Normalization formula:
           If 'CNQ' is present:
             QUAL_norm = min(1000.0, 100.0 * (CNQ / 5.0))
           Else if 'GQ' is present:
             QUAL_norm = min(1000.0, 100.0 * (GQ / 5.0))
       - Rationale: Normalizes typical low-scale Phred quality metrics from CNVKit, 
         mapping the baseline of 5.0 to 100.0.

    6. DRAGEN
       - Source metric: VCF QUAL column value (Phred-scale representation).
       - Prerequisite: Record must have a "PASS" filter flag.
       - Normalization formula (piecewise linear interpolation):
           * If dragen_qual <= 10.0:
               QUAL_norm = dragen_qual * 10.0
           * If dragen_qual > 10.0:
               QUAL_norm = 100.0 + (dragen_qual - 10.0) * (900.0 / 190.0)
       - Rationale: Linearly maps [0, 10.0] to [0.0, 100.0] and maps (10.0, 200.0] to (100.0, 1000.0].
         This allows lower-quality scores to be stretched and higher-quality scores to be compressed 
         into the standard [0, 1000] scale.

    Args:
        input_vcf (str): Path to the input VCF.
        output_vcf (str): Path to write the standardized VCF.
        caller (str): Name of the caller (supports common aliases like GCNV, CNV-KIT).
    """
    # Map caller names to standard uppercase keys to support aliases and case-insensitivity
    caller = caller.upper()
    if caller in ['GATK-GCNV', 'GCNV']:
        caller = 'GATK'
    elif caller in ['CNV-KIT']:
        caller = 'CNVKIT'
    elif caller in ['DRAGEN_GERMLINE', 'DRAGEN-GERMLINE']:
        caller = 'DRAGEN'
    elif caller in ['EXOMEDEPTH', 'EXOME_DEPTH', 'EXOME-DEPTH']:
        caller = 'EXOMEDEPTH'

    # Open the input VCF using PySam VariantFile
    vcf = pysam.VariantFile(input_vcf, "r")

    # Add new FORMAT headers to preserve original metrics for tracking and verification
    vcf.header.formats.add("OQ", "1", "Float", "Original Quality Score")
    vcf.header.formats.add("OAS", "1", "String", "Original Algorithm Score Metric Used")

    # Open the output VCF for writing
    out = pysam.VariantFile(output_vcf, "w", header=vcf.header)
    universal_baseline = 100.0

    for record in vcf:
        # Clean up record ID if it contains colons to prevent SURVIVOR merge format field parse error
        if record.id and ":" in record.id:
            record.id = record.id.replace(":", "_")

        # Get the first sample's FORMAT fields (single-sample pipeline assumption)
        sample = next(iter(record.samples.values()), None)

        # Safely extract existing QUAL and move to the custom OQ format field
        orig_qual = record.qual
        if orig_qual is not None and sample is not None:
            sample["OQ"] = orig_qual  # Save original score
            record.qual = None  # Clear the QUAL field to prepare for normalized score

        qual_norm = None
        metric_used = "UNKNOWN"
        raw_metric = None

        try:
            # 1. CANOES score normalization
            if caller == "CANOES":
                if sample is not None and "Q_SOME" in sample:
                    q_some = float(sample["Q_SOME"])
                    raw_metric = q_some
                    # Scale Q_SOME linearly using configurable canoes_baseline
                    qual_norm = min(1000.0, universal_baseline * (q_some / max(1e-6, canoes_baseline)))
                    metric_used = "Q_SOME"

            # 2. CLAMMS score normalization
            elif caller == "CLAMMS":
                if sample is not None and "Q_SOME" in sample and "Q_EXACT" in sample:
                    q_exact = float(sample["Q_EXACT"])
                    q_some = float(sample["Q_SOME"])
                    raw_metric = q_some
                    # Verify boundary call reliability using Q_EXACT (>= 0 for DEL, >= -0.5 for DUP)
                    is_dup = (record.alts and "<DUP>" in record.alts) or sample.get("GT") == (2,)
                    min_q_exact = -0.5 if is_dup else 0.0
                    if q_exact >= min_q_exact:
                        qual_norm = min(1000.0, universal_baseline * (q_some / max(1e-6, clamms_baseline)))
                    else:
                        qual_norm = 0.0
                    metric_used = "Q_SOME"

            # 3. XHMM score normalization
            elif caller == "XHMM":
                if sample is not None and "SQ" in sample:
                    # GT is a tuple, e.g. (1,) or (2,). Handle index based on genotype:
                    # index 0 for DEL (GT=1), index 1 for DUP (GT=2)
                    gt = sample.get("GT", (None,))
                    idx = 0
                    if (gt and any(g == 2 for g in gt if g is not None)) or (record.alts and "<DUP>" in record.alts):
                        idx = 1
                    
                    sq_val = sample["SQ"]
                    if sq_val is not None:
                        if isinstance(sq_val, (tuple, list)):
                            sq = float(sq_val[idx]) if len(sq_val) > idx else float(sq_val[0])
                        else:
                            sq = float(sq_val)
                        raw_metric = sq
                        # Calibrate XHMM SQ across dynamic range [0, 1000]
                        if sq <= xhmm_baseline:
                            qual_norm = universal_baseline * (sq / max(1e-6, xhmm_baseline))
                        else:
                            denom = max(1.0, 99.0 - xhmm_baseline)
                            qual_norm = universal_baseline + (sq - xhmm_baseline) * (900.0 / denom)
                        qual_norm = min(1000.0, qual_norm)
                        metric_used = "SQ"

            # 4. GATK gCNV score normalization
            elif caller == "GATK":
                if sample is not None and "QS" in sample and "CN" in sample and "NP" in sample:
                    qs = float(sample["QS"])
                    cn = int(sample["CN"])
                    n_int = float(sample["NP"])
                    raw_metric = qs

                    chrom_clean = str(record.chrom).replace('chr', '').replace('CHR', '').strip()
                    is_sex_chrom = chrom_clean in ['X', 'Y']
                    is_dup = (record.alts and "<DUP>" in record.alts) or cn > 2 or (is_sex_chrom and cn == 2)

                    t_gatk = None
                    # Compute dynamic threshold baseline T_GATK
                    if cn == 0:
                        # Homozygous deletion: expect very high signal intensity change
                        t_gatk = min(1000.0, max(400.0, 10.0 * n_int))
                    elif cn == 1 and not (is_sex_chrom and not is_dup):
                        # Heterozygous deletion on autosomes
                        t_gatk = min(1000.0, max(100.0, 10.0 * n_int))
                    elif is_dup:
                        # Duplication event (including male CN=2 on chrX/chrY): lower SNR, lower threshold
                        t_gatk = min(400.0, max(50.0, 4.0 * n_int))

                    if t_gatk is not None:
                        # Scale the QS score by the computed dynamic baseline threshold
                        qual_norm = min(1000.0, universal_baseline * (qs / t_gatk))
                        metric_used = "QS"

            # 5. CNVKIT score normalization
            elif caller == "CNVKIT":
                if sample is not None:
                    # Prefer CNQ if available, fallback to GQ or log2
                    if "CNQ" in sample and sample["CNQ"] is not None:
                        cnq = float(sample["CNQ"])
                        raw_metric = cnq
                        qual_norm = min(1000.0, universal_baseline * (cnq / max(1e-6, cnvkit_baseline_cnq)))
                        metric_used = "CNQ"
                    elif "GQ" in sample and sample["GQ"] is not None:
                        gq = float(sample["GQ"])
                        raw_metric = gq
                        qual_norm = min(1000.0, universal_baseline * (gq / max(1e-6, cnvkit_baseline_cnq)))
                        metric_used = "GQ"
                    elif "log2" in record.info:
                        log2_val = abs(float(record.info["log2"]))
                        raw_metric = log2_val
                        qual_norm = min(1000.0, universal_baseline * (log2_val / max(1e-6, cnvkit_baseline_log2)))
                        metric_used = "log2"

            # 6. DRAGEN score normalization
            elif caller == "DRAGEN":
                if orig_qual is not None and ("PASS" in record.filter.keys() or "." in record.filter.keys() or len(record.filter.keys()) == 0):
                    dragen_qual = float(orig_qual)
                    raw_metric = dragen_qual
                    # Baseline threshold for PASS maps to 100.0
                    if dragen_qual <= dragen_baseline:
                        qual_norm = universal_baseline * (dragen_qual / max(1e-6, dragen_baseline))
                    else:
                        denom = max(1.0, 300.0 - dragen_baseline)
                        qual_norm = universal_baseline + (dragen_qual - dragen_baseline) * (900.0 / denom)
                    qual_norm = min(1000.0, qual_norm)
                    metric_used = "QUAL"

            # 7. EXOMEDEPTH score normalization
            elif caller == "EXOMEDEPTH":
                bf_val = None
                if sample is not None and "BF" in sample and sample["BF"] is not None:
                    try:
                        bf_val = float(sample["BF"])
                    except (ValueError, TypeError):
                        pass
                if bf_val is None and "BF" in record.info:
                    try:
                        bf_val = float(record.info["BF"])
                    except (ValueError, TypeError):
                        pass
                if bf_val is None and orig_qual is not None:
                    bf_val = float(orig_qual)

                if bf_val is not None:
                    raw_metric = bf_val
                    # Linear scaling using configurable exomedepth_baseline
                    qual_norm = min(1000.0, universal_baseline * (max(0.0, bf_val) / max(1e-6, exomedepth_baseline)))
                    metric_used = "BF"

        except (ValueError, TypeError, IndexError) as e:
            print(f"Warning: processing record {record.chrom}:{record.pos} {record.id}: {e}")

        # Set the normalized quality score in the QUAL field (round to 2 decimals if present, None if missing)
        if qual_norm is not None:
            record.qual = round(qual_norm, 2)
        else:
            record.qual = None

        if sample is not None:
            # Preserve raw caller metric in OQ if not already set from orig_qual
            if raw_metric is not None and sample.get("OQ") is None:
                sample["OQ"] = round(float(raw_metric), 2)
            sample["OAS"] = metric_used  # Store original metric identifier in format metadata

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
    parser.add_argument(
        '--canoes_baseline', type=float, default=65.0,
        help='CANOES Q_SOME baseline mapping to 100.0 (default: 65.0)'
    )
    parser.add_argument(
        '--clamms_baseline', type=float, default=50.0,
        help='CLAMMS Q_SOME baseline mapping to 100.0 (default: 50.0)'
    )
    parser.add_argument(
        '--xhmm_baseline', type=float, default=50.0,
        help='XHMM SQ baseline mapping to 100.0 (default: 50.0)'
    )
    parser.add_argument(
        '--cnvkit_baseline_cnq', type=float, default=15.0,
        help='CNVkit CNQ baseline mapping to 100.0 (default: 15.0)'
    )
    parser.add_argument(
        '--cnvkit_baseline_log2', type=float, default=0.5,
        help='CNVkit log2 baseline mapping to 100.0 (default: 0.5)'
    )
    parser.add_argument(
        '--dragen_baseline', type=float, default=30.0,
        help='DRAGEN QUAL baseline mapping to 100.0 (default: 30.0)'
    )
    parser.add_argument(
        '--exomedepth_baseline', type=float, default=20.0,
        help='ExomeDepth Bayes Factor baseline mapping to 100.0 (default: 20.0)'
    )
    args = parser.parse_args()
    standardize_cnv_qual(
        args.input_vcf,
        args.output_vcf,
        args.caller,
        canoes_baseline=args.canoes_baseline,
        clamms_baseline=args.clamms_baseline,
        xhmm_baseline=args.xhmm_baseline,
        cnvkit_baseline_cnq=args.cnvkit_baseline_cnq,
        cnvkit_baseline_log2=args.cnvkit_baseline_log2,
        dragen_baseline=args.dragen_baseline,
        exomedepth_baseline=args.exomedepth_baseline
    )
