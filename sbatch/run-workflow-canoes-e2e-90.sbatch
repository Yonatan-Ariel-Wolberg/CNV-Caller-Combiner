#!/usr/bin/env bash

# =============================================================================
# Detailed documentation: run-workflow-canoes-e2e-90.sbatch
# Purpose: This SLURM submission workflow wrapper supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: End-to-end run path for CANOES.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

# SLURM Job Submission Directives:
# --------------------------------
# --job-name: Distinct identifier for the job in the SLURM queue.
#SBATCH --job-name=combiner-canoes-e2e-90
# --output: Specifies the destination path for standard output (stdout), including job name (%x) and job ID (%j).
#SBATCH --output=slurm-%x-%j.out
# --error: Specifies the destination path for standard error (stderr).
#SBATCH --error=slurm-%x-%j.err
# --partition: The queue partition targeted for execution.
#SBATCH --partition=batch
# --time: Absolute limit on runtime formatted as DD-HH:MM:SS. Set to 72 hours.
#SBATCH --time=72:00:00
# --cpus-per-task: Number of CPU cores allocated per task.
#SBATCH --cpus-per-task=4
# --mem: Total RAM memory limit allocated to the job.
#SBATCH --mem=32G
#SBATCH --exclude=n01,n03

# Shell Configuration:
# -e: Exit immediately if a command exits with a non-zero status.
# -u: Treat unset variables as an error when substituting.
# -o pipefail: The return value of a pipeline is the status of the last command to exit with a non-zero status.
set -euo pipefail

# HPC Module Environment Loading:
# Load required software dependencies using the environment module system.
# - java/21: Required runtime environment for Nextflow.
# - nextflow/25.10: The pipeline orchestrator.
# - samtools/1.23 & bedtools/2.30.0: Sequence alignment utilities.
# - SURVIVOR/1.0.7: SV merging tool.
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7

# Python Environment Activation:
# Activate the python virtual environment containing pandas, matplotlib, and other pipeline dependencies.
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

# Environment Variables:
# Disable ANSI logging features in Nextflow to avoid log file clutter with terminal control characters.
export NXF_ANSI_LOG=false

# Nextflow Pipeline Execution:
# Execute the CANOES end-to-end germline CNV calling workflow.
# Parameters:
# - canoes_germline_cnv_end_to_end.nf: Nextflow script for the CANOES workflow.
# - -profile wits,medium: Config profiles for specific cluster resource settings.
# - --run_calling true: Enables calling process steps.
# - --run_modes "Control,Mosdepth": Specifying control validation and coverage collection mode.
# - --bam_dir: Location of input BAM/CRAM files.
# - --outdir: Output destination for final callsets.
# - --total_samples 90: Specifies the number of target samples in batch.
# - -w: Working directory containing temporary/intermediate task execution folders.
# - -resume: Enables pipeline caching and resumption from the last successfully completed task.
nextflow run canoes_germline_cnv_end_to_end.nf \
  -profile wits,medium \
  --run_calling true \
  --run_modes "Control,Mosdepth" \
  --bam_dir "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/90-samples/stable_bams" \
  --outdir "/home/ywolberg/nf-multi-exome-cnv/results_canoes_e2e_90" \
  --total_samples 90 \
  -w "/home/ywolberg/nf-multi-exome-cnv/work" \
  -resume
