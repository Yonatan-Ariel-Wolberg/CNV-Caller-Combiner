#!/usr/bin/env bash
#SBATCH --job-name=move-results
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

echo "Starting move of results..."
mkdir -p /spaces/ywolberg/cnv-caller-combiner-results/

# Use rsync or mv to move the entire directory safely
rsync -av --remove-source-files /home/ywolberg/cnv-caller-combiner-results/first-test/ /spaces/ywolberg/cnv-caller-combiner-results/first-test/
find /home/ywolberg/cnv-caller-combiner-results/first-test/ -type d -empty -delete || true

echo "Move complete!"
