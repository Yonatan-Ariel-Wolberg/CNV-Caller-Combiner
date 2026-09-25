#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse
import pandas as pd

def run_evaluation(truth_bed, probes_bed, callset_bed):
    """
    Runs evaluate_caller_performance.py on a callset and returns parsed metrics.
    """
    cmd = [
        "python", "bin/evaluate_caller_performance.py",
        "--truth_bed", truth_bed,
        "--probes_bed", probes_bed,
        "--callset_bed", callset_bed
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = {}
        for line in res.stdout.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith("TP:"):
                metrics['TP'] = int(line.split(":")[1].strip())
            elif line.startswith("FP:"):
                metrics['FP'] = int(line.split(":")[1].strip())
            elif line.startswith("FN:"):
                metrics['FN'] = int(line.split(":")[1].strip())
            elif line.startswith("TN:"):
                metrics['TN'] = int(line.split(":")[1].strip())
            elif line.startswith("Sensitivity:"):
                metrics['Sensitivity'] = float(line.split(":")[1].strip())
            elif line.startswith("Precision:"):
                metrics['Precision'] = float(line.split(":")[1].strip())
            elif line.startswith("Specificity:"):
                metrics['Specificity'] = float(line.split(":")[1].strip())
            elif line.startswith("F_beta (beta=2):"):
                metrics['F2_Score'] = float(line.split(":")[1].strip())
            elif line.startswith("Matthews Correlation Coefficient:"):
                metrics['MCC'] = float(line.split(":")[1].strip())
        return metrics
    except Exception as e:
        print(f"Error running evaluation on {callset_bed}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Track CNV caller metrics across pipeline stages.")
    parser.add_argument("--truth_bed", required=True, help="Path to validated truth BED file")
    parser.add_argument("--probes_bed", required=True, help="Path to probes target BED file")
    args = parser.parse_args()

    # Define the 6 progressive stages we want to track
    stages = [
        {"Stage": "1) Raw Callset (Pre-Filter)", "VCF_Pattern": "work/*/*/CANOES_raw_*.vcf"},
        {"Stage": "2) Quality Filtered", "VCF_Pattern": "work/*/*/CANOES_filtered_*.vcf"},
        {"Stage": "3) Intra-Sample Merged", "VCF_Pattern": "work/*/*/merged_intrasample_*.vcf"},
        {"Stage": "4) Inter-Sample Cohort Merged", "VCF_Pattern": "work/*/*/cohort_merged_*.vcf"},
        {"Stage": "5) Rare Variant AF Filtered (<1%)", "VCF_Pattern": "work/*/*/cohort_rare.annotated.tsv"},
        {"Stage": "6) ACMG Stratified (Pathogenic Class 4/5)", "VCF_Pattern": "work/*/*/cohort_acmg_4_5.vcf"}
    ]

    print("Stage Tracking Tool initiated.")
    print("---------------------------------")
    print(f"Truth BED: {args.truth_bed}")
    print(f"Probes BED: {args.probes_bed}\n")

    results = []
    
    # We will simulate/check local patterns
    for s in stages:
        stage_name = s["Stage"]
        pattern = s["VCF_Pattern"]
        matching_files = glob.glob(pattern)
        
        if matching_files:
            # Run evaluation on the first match as representative
            target_file = matching_files[0]
            print(f"Evaluating {stage_name} using: {target_file}")
            
            # Convert VCF to BED if necessary
            temp_bed = "results/temp_stage_eval.bed"
            # Note: in real run, we can parse VCF to BED or call our converter
            # For this tracking tool, we print the file path and placeholder if not converted
            metrics = run_evaluation(args.truth_bed, args.probes_bed, target_file)
            if metrics:
                metrics['Stage'] = stage_name
                results.append(metrics)
        else:
            # Fallback/Placeholder indicating how to query this step
            print(f"No active files found for {stage_name} (Pattern: {pattern}).")

    if results:
        df = pd.DataFrame(results)
        print("\n### PROGRESSION OF METRICS ACROSS STAGES")
        print(df.to_markdown(index=False))

if __name__ == "__main__":
    main()
