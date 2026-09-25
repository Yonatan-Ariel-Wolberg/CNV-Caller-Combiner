#!/bin/bash

# =============================================================================
# Detailed documentation: get_mappability.sh
# Purpose: Downloads and processes UCSC ENCODE Crg mappability BigWig tracks 
#   corresponding to standard read lengths (36, 40, 50, 75, or 100mer) and
#   converts them to CLAMMS-compatible BED format coordinates.
# Author/Context: Part of the CNV-Caller-Combiner pipeline.
# Inputs:
#   - Optional positional argument $1 specifying read length (default: 75).
# Outputs:
#   - mappability.bed: Sorted, CLAMMS-ready, chromosome-stripped BED file.
# Dependencies: bigWigToWig (will download automatically if missing), wget, sed, sort.
# =============================================================================

# Exit immediately if a command exits with a non-zero status, unset variable, or pipe failure
set -euo pipefail

# --- Variables ---
# Read length passed as parameter (defaults to 75 if not provided)
READ_LEN="${1:-75}"
echo "Determined read length: ${READ_LEN}"

# Map read length to the closest UCSC mer (36, 40, 50, 75, 100)
if [ "${READ_LEN}" -le 38 ]; then
    MER=36
elif [ "${READ_LEN}" -le 45 ]; then
    MER=40
elif [ "${READ_LEN}" -le 62 ]; then
    MER=50
elif [ "${READ_LEN}" -le 87 ]; then
    MER=75
else
    MER=100
fi
echo "Selected ${MER}mer mappability track from UCSC."

URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign${MER}mer.bigWig"
BIGWIG_FILE="wgEncodeCrgMapabilityAlign${MER}mer.bigWig"
WIG_FILE="wgEncodeCrgMapabilityAlign${MER}mer.wig"
OUTPUT_BED="mappability.bed"

# --- Step 1: Check Dependencies ---
echo "Checking for required tools..."
if ! command -v bigWigToWig &> /dev/null; then
    if [ -f "./bigWigToWig" ]; then
        echo "Found bigWigToWig in current directory."
        export PATH=".:$PATH"
    else
        echo "'bigWigToWig' is not installed or not in your PATH. Downloading from UCSC..."
        wget -qO bigWigToWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
        chmod +x bigWigToWig
        export PATH=".:$PATH"
    fi
fi

# --- Step 2: Download ---
echo "Downloading mappability track (this may take a while)..."
# Only download if the file doesn't already exist
if [ ! -f "$BIGWIG_FILE" ]; then
    wget "$URL"
else
    echo "Found existing $BIGWIG_FILE, skipping download."
fi

# --- Step 3: Convert to WIG ---
echo "Converting bigWig to WIG format..."
bigWigToWig "$BIGWIG_FILE" "$WIG_FILE"

# --- Step 4: Parse to BED ---
echo "Formatting to CLAMMS-ready BED file..."
# grep -v '^#' removes header lines
# sed 's/^chr//g' removes the 'chr' prefix from chromosome names
# LC_ALL=C sort -k1,1 -k2,2n sorts the BED file deterministically by chromosome and start position
grep -v '^#' "$WIG_FILE" | sed 's/^chr//g' | LC_ALL=C sort -k1,1 -k2,2n > "$OUTPUT_BED"

echo "Success! Your CLAMMS-ready file is saved as: $OUTPUT_BED"

# Optional: Clean up large intermediate files to save disk space
# Remove the '#' from the line below to enable automatic cleanup
# rm "$BIGWIG_FILE" "$WIG_FILE"
