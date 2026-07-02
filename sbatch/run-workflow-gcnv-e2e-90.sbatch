#!/usr/bin/env bash

# =============================================================================
# Detailed documentation: run-workflow-gcnv-e2e-90.sbatch
# Purpose: This SLURM submission workflow wrapper supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: End-to-end run path for GATK-gCNV.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
#SBATCH --job-name=combiner-gcnv-e2e-90
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# Load environment modules (HPC Native)
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7

# Activate python virtual environment
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

export NXF_ANSI_LOG=false

mkdir -p /home/ywolberg/run_gcnv
cd /home/ywolberg/run_gcnv
cp -r /home/ywolberg/nf-multi-exome-cnv/.nextflow .
find .nextflow/ -name "LOCK" -delete 2>/dev/null || true

nextflow run /home/ywolberg/CNV-Caller-Combiner/gcnv_germline_cnv_end_to_end.nf \
  -profile wits,medium \
  --run_calling true \
  --run_modes "Control,Mosdepth" \
  --bam_dir "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/90-samples/stable_bams" \
  --outdir "/home/ywolberg/nf-multi-exome-cnv/results_gcnv_e2e_90" \
  --total_samples 90 \
  -w "/home/ywolberg/nf-multi-exome-cnv/work" \
  -resume
