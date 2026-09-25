#!/usr/bin/env bash
#SBATCH --job-name=run-hard-reruns-combined
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --exclude=n01,n03

set -euo pipefail

# HPC Module Environment Loading:
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7 || true

# Execute the master rerun script
bash /home/ywolberg/CNV-Caller-Combiner/scratch/run_hard_reruns_combined.sh
