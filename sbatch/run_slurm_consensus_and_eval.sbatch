#!/usr/bin/env bash
#SBATCH --job-name=consensus-and-eval
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# Load environment modules
module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7

# Python Environment Activation (needed for evaluation scripts inside e2e)
source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate

export NXF_ANSI_LOG=false

# Clean stale locks in the work directory if any
find /home/ywolberg/cnv-caller-combiner-results/first-test/work/ -name "LOCK" -delete 2>/dev/null || true

cd /home/ywolberg/CNV-Caller-Combiner

# =====================================================================
# 1. Run the consensus workflows (with AnnotSV 3.5.10 auto-download)
# =====================================================================
echo "Running SURVIVOR Control consensus..."
nextflow run cnv-caller-combiner.nf \
  -profile medium \
  --workflow survivor \
  --canoes_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/canoes \
  --clamms_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/clamms \
  --xhmm_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/xhmm \
  --cnvkit_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/cnvkit \
  --gcnv_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/gcnv \
  --outdir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_survivor_control \
  -resume

echo "Running Truvari Control consensus..."
nextflow run cnv-caller-combiner.nf \
  -profile medium \
  --workflow truvari \
  --canoes_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/canoes \
  --clamms_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/clamms \
  --xhmm_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/xhmm \
  --cnvkit_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/cnvkit \
  --gcnv_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Control/gcnv \
  --outdir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_truvari_control \
  -resume

echo "Running SURVIVOR Mosdepth consensus..."
nextflow run cnv-caller-combiner.nf \
  -profile medium \
  --workflow survivor \
  --canoes_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/canoes \
  --clamms_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/clamms \
  --xhmm_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/xhmm \
  --cnvkit_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/cnvkit \
  --gcnv_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/gcnv \
  --outdir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_survivor_mosdepth \
  -resume

echo "Running Truvari Mosdepth consensus..."
nextflow run cnv-caller-combiner.nf \
  -profile medium \
  --workflow truvari \
  --canoes_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/canoes \
  --clamms_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/clamms \
  --xhmm_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/xhmm \
  --cnvkit_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/cnvkit \
  --gcnv_dir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_inputs/Mosdepth/gcnv \
  --outdir /home/ywolberg/cnv-caller-combiner-results/first-test/CNV-Caller-Combiner-results/consensus_truvari_mosdepth \
  -resume

# =====================================================================
# 2. Run the comprehensive single-caller evaluations
# =====================================================================
echo "Running comprehensive single-caller evaluations..."
python3 /home/ywolberg/.gemini/antigravity-cli/brain/90368432-f5c4-4ca0-9b62-3c02b22be24c/scratch/evaluate_all_callers_comprehensive.py

echo "SLURM Batch Job Completed Successfully!"
