#!/usr/bin/env bash
#SBATCH --job-name=run-all-tests-resume
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --exclude=n01,n03

set -euo pipefail

# Output and work directories in symlinked path
mkdir -p /home/ywolberg/cnv-caller-combiner-results/first-test/work
mkdir -p /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results
mkdir -p /home/ywolberg/cnv-caller-combiner-results/first-test/nf-exomecnv-results

# Load modules
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7

# Python Environment Activation (needed for evaluation scripts inside e2e)
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

export NXF_ANSI_LOG=false

# Clean stale locks in the work directory if any
find /home/ywolberg/cnv-caller-combiner-results/first-test/work/ -name "LOCK" -delete 2>/dev/null || true

# =========================================
# 1. Run e2e workflows in CNV-Caller-Combiner
# =========================================
cd /home/ywolberg/CNV-Caller-Combiner

for tool in canoes clamms cnvkit gcnv xhmm; do
  echo "========================================="
  echo "Running e2e workflow for ${tool}..."
  echo "========================================="
  nextflow run /home/ywolberg/CNV-Caller-Combiner/${tool}_germline_cnv_end_to_end.nf \
    -profile wits,medium \
    --run_calling true \
    --run_modes "Control,Mosdepth" \
    --bam_dir "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/90-samples/stable_bams" \
    --outdir "/home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/${tool}" \
    --total_samples 90 \
    -w "/home/ywolberg/cnv-caller-combiner-results/first-test/work" \
    -resume
done

# =========================================
# 2. Run subworkflows in nf-exomecnv
# =========================================
cd /home/ywolberg/nf-exomecnv

for tool in canoes clamms xhmm; do
  echo "========================================="
  echo "Running nf-exomecnv subworkflow for ${tool}..."
  echo "========================================="
  nextflow run /home/ywolberg/nf-exomecnv/main.nf \
    -profile wits \
    -params-file /home/ywolberg/nf-exomecnv/params-${tool}-wits-ddd-uk-90.json \
    --outdir "/home/ywolberg/cnv-caller-combiner-results/first-test/nf-exomecnv-results/${tool}" \
    -w "/home/ywolberg/cnv-caller-combiner-results/first-test/work" \
    -resume
done

echo "All tasks finished!"
