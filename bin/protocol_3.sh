#!/bin/bash

# =============================================================================
# Detailed documentation: protocol_3.sh
# Purpose: This pipeline script identifies rare and putative de novo CNV calls using
#   a combination of PLINK (v1.07/v1.9) and PLINK/Seq tools. It filters out common
#   CNVs (frequency > 10%), loads the remaining variants and a pedigree file into
#   a PLINK/Seq project, and runs quality-based filtering to detect de novo events.
#
# Arguments:
#   - $1: Input VCF file containing candidate CNV calls.
#
# Environment Variables & External Dependencies:
#   - xcnv_to_cnv: Utility to convert XHMM .xcnv output format to PLINK .cnv format.
#   - /home/phelelani/applications/plink/plink: PLINK executable.
#   - /home/phelelani/applications/plinkseq/pseq: PLINK/Seq executable.
#   - DATA.xcnv: Pre-existing XHMM CNV calls.
#   - DATA.pedinfo: Pre-existing pedigree information file.
#   - EXOME.interval_list: Pre-existing exome bait/target intervals.
#
# Outputs:
#   - DATA.fam: PLINK family description file constructed from the VCF sample headers.
#   - DATA.cnv: XHMM calls formatted as PLINK CNVs.
#   - DATA.maf_0.1.vcf: VCF file filtered to exclude common CNV regions.
#   - DATA.maf_0.1.CNV_regions.txt: List of coordinates for excluded common CNV regions.
#   - PLINK/Seq project structure (`new-project` containing variant database).
#   - annotated_targets.refseq: Intersected targets with RefSeq annotations.
# =============================================================================

# Argument: Input VCF file containing raw CNV calls.
VCF=$1

## 1. FILTER CNV BASED ON FREQUENCY ESTIMATED FROM THE SAMPLE
##-- SCRIPT TO CONVERT XHMM CNV CALLS TO PLINK FORMAT (REQUIRES PLINK 1.7)
# Extract sample names from VCF header (columns 10 to end) and write to a PLINK .fam file.
# Each sample is assigned family ID, individual ID, paternal ID (0), maternal ID (0), sex (1=male), and phenotype (1=unaffected).
grep "^#CHROM" $VCF | awk '{for (i = 10; i <= NF; i++) print $i,1,0,0,1,1}' > DATA.fam

# Convert XHMM raw CNV calls (DATA.xcnv) into PLINK CNV format (DATA.cnv).
xcnv_to_cnv DATA.xcnv > DATA.cnv

# Generate a map of CNV regions across samples using PLINK.
/home/phelelani/applications/plink/plink --noweb --cfile DATA --cnv-make-map --out DATA

##-- FIND ANY CNV THAT OVERLAPS 50% OF ANOTHER CNV THAT OCCURS IN MORE THAN 10% OF ALL SAMPLES
##-- (A RELATIVELY LIBERAL FREQUENCY THRESHOLD TO REMOVE ONLY COMMON CNV AND ARTIFACTS).
## CALCULATE THE NUMBER OF SAMPLES IN THE DATA
NUM_SAMPLES=`cat DATA.fam | wc -l`

## SET THE FREQUENCY THRESHOLD AT 10% OF THE TOTAL NUMBER OF SAMPLES
# Calculate 10% ceiling threshold of samples to identify common CNVs.
THRESH_NUM=`echo "0.1 * $NUM_SAMPLES" | bc | awk '{x = $1; if (x != int(x)) {x = int(x)+1} print x}'`

# Exclude CNVs occurring in more than 10% of samples and overlap by 50% or more.
/home/phelelani/applications/plink/plink --noweb --cfile DATA --cnv-freq-exclude-above $THRESH_NUM --cnv-overlap 0.5 --cnv-write --out DATA.maf_0.1

# Format common CNV regions as 'chr:start-end' coordinates.
cat DATA.maf_0.1.cnv | awk '{if (NR>1) print $3":"$4"-"$5}' | sort | uniq | sed 's/^/chr/g' > DATA.maf_0.1.CNV_regions.txt

## EXCLUDE ALL SUCH COMMON REGIONS FROM THE VCF FILE
# Loop through the VCF records and omit ones that match an ID/Region in the common CNV list.
cat $VCF | awk '{if (substr($1,1,1)=="#") {print $_} else {found=-1; cmd = "grep -w "$3" DATA.maf_0.1.CNV_regions.txt"; cmd | getline found; close(cmd);
if (found != -1) {print $_}}}' > DATA.maf_0.1.vcf

## 2. DETECT PUTATIVE DE NOVO CNV BY STRICT USE OF XHMM QUALITY SCORES
## CREATE A NEW PLINK/SEQ PROJECT THAT WILL CONTAIN THE XHMM-CALLED VARIANTS:
/home/phelelani/applications/plinkseq/pseq DATA new-project

## LOAD A PLINK/SEQ PEDINFO FILE CONTAINING THE CHILD-PARENT PEDIGREE RELATIONSHIPS (NAMED ‘DATA.PEDINFO’ HERE):
/home/phelelani/applications/plinkseq/pseq DATA load-pedigree --file DATA.pedinfo

## LOAD THE FREQUENCY-FILTERED VCF FILE INTO THE PLINK/SEQ PROJECT:
/home/phelelani/applications/plinkseq/pseq DATA load-vcf --vcf DATA.maf_0.1.vcf

## NOW, ASSUMING YOU'VE CHOSEN A QUALITY SCORE FILTER OF 60 TO CALL DE NOVO CNV ON THE AUTOSOMES, RUN:
# Filters out sex chromosomes and requires both site quality (SQ) and individual genotype quality (NQ) >= 60.
/home/phelelani/applications/plinkseq/pseq DATA cnv-denovo --mask reg.ex=chrX,chrY --minSQ 60 --minNQ 60 --out DATA

## PROTOCOL 2 (Intersecting targets with RefSeq annotations)
# Intersect exome target intervals with RefSeq groups to map target exons to annotated genes.
pseq . loc-intersect --group refseq --locdb locdb --file EXOME.interval_list --out annotated_targets.refseq
