#!/bin/bash

# =============================================================================
# Detailed documentation: combine_picard_qc_metrics.sh
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

## GET FILE HEADERS
grep --no-filename --no-group-separator -A1 "## METRICS CLASS" *.hs_metrics.txt | grep -v "## METRICS CLASS" | uniq | awk -v OFS="\t" '{ print "SAMPLE",$9,$32,$48,$52,$63,$64 }' > hs_metrics.txt
grep --no-filename --no-group-separator -A1 "## METRICS CLASS" *.insert_size_metrics.txt | grep -v "## METRICS CLASS" | uniq | awk  -v OFS="\t" '{ print "SAMPLE",$6}' > insert_size_metrics.txt

## PRINT OUT FOR HS METRICS
for hs_metrics in *.hs_metrics.txt
do
    sample="${hs_metrics%.hs_metrics.txt}"
    sed -n '/^BAIT_SET\(.*\)$/,/^$/p' $hs_metrics |\
        sed '1d;$d' | awk -v OFS="\t" -v sample="$sample" '{ print sample,$9,$32,$48,$52,$63,$64 }'
done >> hs_metrics.txt

## PRINT OUT FOR INSERT SIZE METRICS
for ins_metrics in *.insert_size_metrics.txt
do
    sample="${ins_metrics%.insert_size_metrics.txt}"
    sed -n '/^MEDIAN_INSERT_SIZE\(.*\)$/,/^$/p' $ins_metrics |\
        sed '1d;$d' | awk -v OFS="\t" -v sample="$sample" '{ print sample,$6}'
done >> insert_size_metrics.txt

awk -F'\t' -v OFS='\t' 'NR==FNR {a[$1]=$0; next} $1 in a {print a[$1], $2}' hs_metrics.txt insert_size_metrics.txt > qcs_metrics
