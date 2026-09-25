#!/usr/bin/env bash

# =============================================================================
# SBATCH SCRIPT: RERUN SUITE A5 & A6 UNDER UPDATED HIGH PRECISION PROFILE
# =============================================================================
# Purpose:
#   Re-executes Suite A5 (Control - High Precision Truvari) and Suite A6
#   (Control - High Precision SURVIVOR) using the updated distinct and
#   strictly monotonic CLAMMS threshold (clamms_q_exact_dup_min = -0.2).
#   Following completion, recalculates the global comparative benchmark
#   statistical report across all 24 pipeline branches.
# =============================================================================

#SBATCH --job-name=rerun-a5-a6
#SBATCH --output=/spaces/ywolberg/qc-gated-cnv-combiner/slurm-rerun-a5-a6-%j.out
#SBATCH --error=/spaces/ywolberg/qc-gated-cnv-combiner/slurm-rerun-a5-a6-%j.err
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --exclude=n18,n26,n28

set -euo pipefail

BASE_OUTDIR="/spaces/ywolberg/qc-gated-cnv-combiner"
WORKDIR="/spaces/ywolberg/2-test-cnv-caller-combiner/work"
CACHE_DIR="/spaces/ywolberg/nxf_cache"
QC_PASSED_SAMPLESHEET="${BASE_OUTDIR}/samplesheet_bams_qc_passed.tsv"

export NXF_SINGULARITY_CACHEDIR="${CACHE_DIR}"
export NXF_ANSI_LOG=false

module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7 || true
if [ -f "/home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate" ]; then
    source /home/ywolberg/dragen-germline-enrichment-v4-4-6/ddd-test/100x-ddduk-run/filter-merge-freq/.venv/bin/activate
fi

cd /home/ywolberg/CNV-Caller-Combiner

# Ensure converted_bams symlinks are intact in case storage cleanup ran
mkdir -p "${BASE_OUTDIR}/converted_bams"
for bam in /spaces/ywolberg/2-test-cnv-caller-combiner/converted_bams/*.bam; do
    [ -f "$bam" ] || continue
    fname=$(basename "$bam")
    [ -e "${BASE_OUTDIR}/converted_bams/${fname}" ] || ln -sf "$bam" "${BASE_OUTDIR}/converted_bams/${fname}"
    [ -e "${BASE_OUTDIR}/converted_bams/${fname}.bai" ] || ln -sf "${bam}.bai" "${BASE_OUTDIR}/converted_bams/${fname}.bai"
done

run_rerun_phase() {
    local phase_name="$1"
    local params_file="$2"
    local merger="$3"
    local out_subdir="$4"
    local profile="$5"

    echo "=========================================================================="
    echo " RUNNING: ${phase_name}"
    echo " Merger: ${merger} | Profile: ${profile}"
    echo " Start: $(date)"
    echo "=========================================================================="
    python3 /home/ywolberg/CNV-Caller-Combiner/bin/nxf_lock_manager.py clean || true

    nextflow run cnv-caller-combiner.nf \
      -profile "wits,medium,${profile}" \
      -params-file "${params_file}" \
      --workflow local \
      --merger_mode "${merger}" \
      --samplesheet_bams "${QC_PASSED_SAMPLESHEET}" \
      --outdir "${BASE_OUTDIR}/${out_subdir}" \
      -w "${WORKDIR}" \
      -resume
}

echo "Starting High Precision Re-Execution for Suite A (Arms 9, 10, 11, 12)..."

# --- Suite A5 Rerun: Control - High Precision (Truvari) ---
run_rerun_phase "Suite A5 (Rerun): Control - High Precision (Truvari)" \
  "params/ddd-uk/params-local-2-test-control-high-precision.json" "truvari" "control/high_prec_truvari" "high_precision"

# --- Suite A6 Rerun: Control - High Precision (SURVIVOR) ---
run_rerun_phase "Suite A6 (Rerun): Control - High Precision (SURVIVOR)" \
  "params/ddd-uk/params-local-2-test-control-high-precision.json" "survivor" "control/high_prec_survivor" "high_precision"

echo "=========================================================================="
echo " Recalculating Comparative Statistical Significance Across All 24 Arms"
echo " Start Time: $(date)"
echo "=========================================================================="

mkdir -p "${BASE_OUTDIR}/benchmarks"
python3 /home/ywolberg/CNV-Caller-Combiner/bin/compare_pipeline_benchmarks.py \
    --input-eval-dir "${BASE_OUTDIR}" \
    --output-tsv "${BASE_OUTDIR}/benchmarks/pipeline_statistical_comparisons.tsv" \
    --output-report "${BASE_OUTDIR}/benchmarks/pipeline_statistical_significance_report.md" \
    --bootstrap-iterations 1000 \
    --seed 42 || echo "Warning: Statistical comparison engine encountered an issue."

echo "=========================================================================="
echo " High Precision Re-Execution & Statistical Harmonization Finished!"
echo " End Time: $(date)"
echo "=========================================================================="
