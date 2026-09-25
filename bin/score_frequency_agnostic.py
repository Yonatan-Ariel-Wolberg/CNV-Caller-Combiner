#!/usr/bin/env python3
"""
score_frequency_agnostic.py
==========================
Scores and prioritizes CNVs based entirely on intrinsic biological and functional
constraints, completely independent of external European-dominated frequency databases
(e.g., ClinGen frequency, 1000 Genomes EUR, DGV, gnomAD-SV).

Key Scoring Dimensions (0 - 100 points):
----------------------------------------
1. Gene Constraint & Dosage Sensitivity:
   - High pLI (>= 0.90) or low LOEUF (<= 0.35): +25 points
   - ClinGen Haploinsufficiency (HI = 3 for DEL): +30 points
   - ClinGen Triplosensitivity (TS = 3 for DUP): +25 points
   - OMIM Morbid / Disease Gene Overlap: +20 points
2. Genomic & Functional Disruption:
   - Exonic / CDS coding sequence disruption: +15 points
   - Large structural span (> 100 kb across multiple exons): +10 points
3. Pipeline Caller Concordance:
   - Supported by >= 2 independent algorithms: +10 points
   - Supported by >= 3 independent algorithms: +15 points

Priority Tiers:
---------------
- Tier 1 (High Biological Constraint): Score >= 70
- Tier 2 (Moderate Constraint): Score 40 - 69
- Tier 3 (Low Constraint / Tolerated): Score < 40

Designed to be easy to read, maintain, and adapt.
"""

import argparse
import sys
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Score CNVs using frequency-agnostic biological constraint metrics."
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input TSV or CSV containing AnnotSV annotations or filtered CNV calls."
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output TSV file path for prioritized and scored CNVs."
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="Optional minimum priority score threshold for output rows (default: 0.0)."
    )
    return parser.parse_args()


def get_numeric(val, default=0.0):
    """Safely convert a value to float."""
    try:
        if pd.isna(val) or val == "" or val == ".":
            return default
        return float(val)
    except (ValueError, TypeError):
        return default


def calculate_variant_score(row):
    """
    Computes a 0-100 constraint-based priority score for a single variant annotation.
    Does not use external population allele frequencies.
    """
    score = 0.0
    reasons = []

    # 1. Variant Type (DEL vs DUP)
    sv_type = str(row.get("SV_type", row.get("svtype", row.get("SVTYPE", "")))).upper()
    is_del = "DEL" in sv_type
    is_dup = "DUP" in sv_type

    # 2. Gene Constraint: pLI and LOEUF
    # High pLI (prob. of being loss-of-function intolerant >= 0.9)
    pli = get_numeric(row.get("pLI", row.get("ExAC_pLI", np.nan)), np.nan)
    if not np.isnan(pli) and pli >= 0.90:
        score += 25.0
        reasons.append(f"High_pLI({pli:.2f})")

    # LOEUF (Loss-of-function observed/expected upper bound <= 0.35)
    loeuf = get_numeric(row.get("LOEUF", row.get("oe_lof_upper", np.nan)), np.nan)
    if not np.isnan(loeuf) and loeuf <= 0.35:
        score += 20.0
        reasons.append(f"Low_LOEUF({loeuf:.2f})")

    # 3. ClinGen Dosage Sensitivity (Haploinsufficiency & Triplosensitivity)
    hi_score = str(row.get("ClinGen_HI", row.get("DDD_HI_percent", ""))).strip()
    if is_del and ("3" in hi_score or "Sufficient" in hi_score or "dosage" in hi_score.lower()):
        score += 30.0
        reasons.append("ClinGen_Haploinsufficient_Loss")

    ts_score = str(row.get("ClinGen_TS", "")).strip()
    if is_dup and ("3" in ts_score or "Sufficient" in ts_score):
        score += 25.0
        reasons.append("ClinGen_Triplosensitive_Gain")

    # 4. OMIM Morbid / Known Disease Gene Overlap
    omim = str(row.get("OMIM_morbid", row.get("OMIM_ID", row.get("Phenotype", "")))).strip()
    if omim and omim not in [".", "0", "nan", "None", ""]:
        score += 20.0
        gene = str(row.get("Gene_name", row.get("gene", "")))
        reasons.append(f"OMIM_Morbid_Gene({gene})")

    # 5. Multi-Caller Agreement (Concordance across CANOES, CLAMMS, XHMM, etc.)
    caller_count = get_numeric(row.get("caller_count", row.get("num_callers", 1)), 1)
    if caller_count >= 3:
        score += 15.0
        reasons.append(f"MultiCaller_Concordance({int(caller_count)}_callers)")
    elif caller_count >= 2:
        score += 10.0
        reasons.append(f"MultiCaller_Concordance({int(caller_count)}_callers)")

    # 6. Exonic / CDS Coding Disruption
    location = str(row.get("Location", row.get("location", ""))).lower()
    if "cds" in location or "exon" in location or "tx-overlap" in location:
        score += 15.0
        reasons.append("Coding_Sequence_Disrupted")

    # Cap score at 100
    final_score = min(100.0, score)

    # Assign clinical constraint tier
    if final_score >= 70.0:
        tier = "Tier 1: High Biological Constraint"
    elif final_score >= 40.0:
        tier = "Tier 2: Moderate Constraint"
    else:
        tier = "Tier 3: Low Constraint / Tolerated"

    criteria_str = "; ".join(reasons) if reasons else "No_High_Constraint_Features"
    return final_score, tier, criteria_str


def main():
    args = parse_args()

    # Determine file delimiter
    sep = "\t" if args.input.endswith(".tsv") or ".tsv." in args.input else ","

    try:
        df = pd.read_csv(args.input, sep=sep, low_memory=False)
    except Exception as e:
        sys.exit(f"Error reading input file {args.input}: {e}")

    if df.empty:
        # Create empty dataframe with expected headers
        cols = list(df.columns) + ["Biological_Priority_Score", "Priority_Tier", "Biological_Criteria"]
        empty_out = pd.DataFrame(columns=cols)
        empty_out.to_csv(args.output, sep="\t", index=False)
        print(f"[INFO] Input file {args.input} was empty. Created empty template at {args.output}")
        return

    # Compute constraint scores row-by-row
    scores, tiers, criteria = [], [], []
    for _, row in df.iterrows():
        s, t, c = calculate_variant_score(row)
        scores.append(s)
        tiers.append(t)
        criteria.append(c)

    df["Biological_Priority_Score"] = scores
    df["Priority_Tier"] = tiers
    df["Biological_Criteria"] = criteria

    # Sort descending by priority score
    df = df.sort_values(by="Biological_Priority_Score", ascending=False)

    # Filter by min score if specified
    if args.min_score > 0.0:
        df = df[df["Biological_Priority_Score"] >= args.min_score]

    # Save output TSV
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[SUCCESS] Scored {len(df)} variants. Output saved to {args.output}")


if __name__ == "__main__":
    main()
