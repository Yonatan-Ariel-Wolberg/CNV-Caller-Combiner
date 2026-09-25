#!/usr/bin/env bash

#SBATCH --job-name=combiner-clamms-e2e-90-hard
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --exclude=n01,n03

set -euo pipefail

# HPC Module Environment Loading:
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7

# Python Environment Activation:
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

# Environment Variables:
export NXF_ANSI_LOG=false

# Sandboxed Run Directory Preparation:
mkdir -p /home/ywolberg/run_clamms
cd /home/ywolberg/run_clamms
cp -r /home/ywolberg/nf-multi-exome-cnv/.nextflow . 2>/dev/null || true
find .nextflow/ -name "LOCK" -delete 2>/dev/null || true

# Nextflow Pipeline Execution:
nextflow run /home/ywolberg/CNV-Caller-Combiner/clamms_germline_cnv_end_to_end.nf \
  -profile wits,medium \
  --run_calling true \
  --run_modes "Control,Mosdepth" \
  --bam_dir "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/90-samples/stable_bams" \
  --outdir "/home/ywolberg/nf-multi-exome-cnv/results_clamms_e2e_90" \
  --total_samples 90 \
  -w "/home/ywolberg/nf-multi-exome-cnv/work" \
  -resume
