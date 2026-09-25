#!/bin/bash
#SBATCH --job-name=gcnv_quick_test
#SBATCH --output=run_quick_test.out
#SBATCH --error=run_quick_test.err
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Define paths
PROBES_CHR20_22="/home/ywolberg/CNV-Caller-Combiner/scratch/probes_chr20_22.bed"
PARAMS_FILE="/home/ywolberg/CNV-Caller-Combiner/params/ddd-uk/params-local-quick-test.json"

echo "Running quick test with GATK gCNV model training enabled (Chromosomes 20-22)"
echo "Combination: QC=true, BED=chr20-22, Merger=Truvari"

nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc true \
  --probes_bed $PROBES_CHR20_22 \
  --merger_mode truvari \
  -resume
