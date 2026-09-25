#!/usr/bin/env bash

# =============================================================================
# SBATCH Script to run quality assessment on 90 DDD-UK samples
# =============================================================================

#SBATCH --job-name=qc-90
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=n01,n03

set -euo pipefail

# Load required cluster modules
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

export NXF_ANSI_LOG=false

# Create and navigate to a separate launch directory to avoid session lock conflicts
mkdir -p /home/ywolberg/CNV-Caller-Combiner/qc_launch
cd /home/ywolberg/CNV-Caller-Combiner/qc_launch

# Run target-aware quality assessment workflow
nextflow run /home/ywolberg/CNV-Caller-Combiner/target-aware-exome-sequencing-quality-assess.nf \
  -c /home/ywolberg/CNV-Caller-Combiner/nextflow.config \
  -profile wits,medium \
  --input_dir "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/90-samples/stable_bams" \
  --ref_genome "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --target_bed "/dataG/ddd/data/resources/canoes/probes_sanger.bed" \
  --output_dir "/home/ywolberg/CNV-Caller-Combiner/qc_results_ddd_uk" \
  -resume
