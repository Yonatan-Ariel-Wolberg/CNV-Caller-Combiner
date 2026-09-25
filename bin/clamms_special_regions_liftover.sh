#!/bin/bash

# =============================================================================
# Detailed documentation: clamms_special_regions_liftover.sh
# Purpose: This script downloads the original hg19 special regions file from CLAMMS,
#   downloads UCSC liftOver tool, downloads hg19 to hg38 chain file, converts the
#   coordinates to hg38, and sorts the output BED file for compatibility.
# Author/Context: Part of the CNV-Caller-Combiner pipeline.
# Inputs: None.
# Outputs:
#   - clamms_special_regions_hg38.bed: Mapped/lifted over coordinates to hg38.
#   - unmapped_regions.bed: Coordinate regions failed to liftover.
# Dependencies: wget, chmod, sort.
# =============================================================================

# Exit immediately if a command fails, unset variable, or pipe failure
set -euo pipefail

# --- Step 1: Download the original hg19 file from CLAMMS ---
echo "Downloading clamms_special_regions.bed (hg19)..."
wget -qO clamms_special_regions_hg19.bed https://raw.githubusercontent.com/rgcgithub/clamms/master/data/clamms_special_regions.bed

# --- Step 2: Download the UCSC liftOver tool ---
echo "Downloading liftOver tool..."
wget -qO liftOver http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

# --- Step 3: Download the hg19 to hg38 chain file ---
echo "Downloading chain file..."
wget -qO hg19ToHg38.over.chain.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# --- Step 4: Run the conversion ---
echo "Converting coordinates to hg38..."
./liftOver \
  clamms_special_regions_hg19.bed \
  hg19ToHg38.over.chain.gz \
  clamms_special_regions_hg38.tmp.bed \
  unmapped_regions.bed

# Sort the liftover output by chromosome and numeric start coordinate
LC_ALL=C sort -k1,1 -k2,2n clamms_special_regions_hg38.tmp.bed > clamms_special_regions_hg38.bed
rm -f clamms_special_regions_hg38.tmp.bed

echo "Success! Your hg38 file is saved as: clamms_special_regions_hg38.bed"
