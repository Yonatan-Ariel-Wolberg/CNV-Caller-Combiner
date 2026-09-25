#!/bin/bash

# =============================================================================
# Detailed documentation: custom_ref_panel.sh
# Purpose: This pipeline helper script builds a custom reference panel for CLAMMS
#   (Copy Number Estimation using Multi-Sample Sequencing data). It aggregates
#   normalized coverage data across all sample-level `.norm.cov.bed` files, structures
#   them into a matrix, runs Singular Value Decomposition (SVD) to perform PCA,
#   and outputs PCA coordinates for nearest-neighbor cohort selection.
#
# Environment Variables Required:
#   - $CLAMMS_DIR: Path to the CLAMMS installation directory containing transpose.gawk.
#
# External Command Dependencies:
#   - svd (from SVDLIBC or CLAMMS package): Performs singular value decomposition.
#   - gawk: GNU awk used for running the transpose scripts.
#
# Inputs:
#   - *.norm.cov.bed: Normalized coverage BED files for all samples in the directory.
#
# Outputs:
#   - matrix.txt: Intermediate file with matrix header and coverage values.
#   - svd-output-*: SVD execution outputs (specifically svd-output-Ut).
#   - sample.names.txt: List of sample names.
#   - pca.coordinates.txt: Tab-separated coordinates of the sample PCA mapping.
# =============================================================================

# Compute the total number of normalized coverage BED files (representing samples).
NUM_SAMPLES=`ls *.norm.cov.bed | wc -l | awk '{print $1}'`

# Compute the number of autosome windows (genomic intervals) across BED files.
# Exclude sex chromosomes (X/Y) and mitochondrial contigs (M/MT).
FIRST_BED=$(ls *.norm.cov.bed | head -n 1)
NUM_WINDOWS=$(awk '$1 != "X" && $1 != "Y" && $1 != "chrX" && $1 != "chrY" && $1 != "M" && $1 != "chrM" && $1 != "MT" && $1 != "chrMT" {x++} END {print x}' "${FIRST_BED}")

# Create the matrix.txt header containing the dimensions: (NUM_SAMPLES x NUM_WINDOWS)
echo -e "$NUM_SAMPLES\t$NUM_WINDOWS" > matrix.txt

# Iterate through each normalized coverage file.
# Extract the normalized coverage values (4th column) for all autosome intervals,
# transpose them to a single row using CLAMMS transpose utility, and append to matrix.txt.
ls *.norm.cov.bed | while read FILE
do
    awk '$1 != "X" && $1 != "Y" && $1 != "chrX" && $1 != "chrY" && $1 != "M" && $1 != "chrM" && $1 != "MT" && $1 != "chrMT" { print $4 }' "${FILE}" | gawk -f $CLAMMS_DIR/transpose.gawk >> matrix.txt
done

# Run SVD with 4 dimensions on matrix.txt to extract eigenvectors (e.g. principal components).
# Options:
#   -d 4: compute 4 singular triplets
#   -o svd-output: output prefix
#   -r dt: input is in text format
svd -d 4 -o svd-output -r dt matrix.txt

# Extract sample names by stripping the `.norm.cov.bed` extension.
ls *.norm.cov.bed | cut -d '.' -f 1 > sample.names.txt

# Format SVD output:
#   - tail -n +2: skip the header of the Ut matrix.
#   - tr ' ' '\t': convert space delimiters to tabs.
#   - gawk -f ...: transpose the Ut matrix to align dimensions.
#   - paste: combine the sample names column with the PCA coordinates column.
tail -n +2 svd-output-Ut | tr ' ' '\t' | gawk -f $CLAMMS_DIR/transpose.gawk | paste sample.names.txt - > pca.coordinates.txt
