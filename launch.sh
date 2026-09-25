#!/usr/bin/env bash
# ==============================================================================
# launch.sh
#
# Production & Beginner-Friendly Unified Launcher for CNV-Caller-Combiner.
# Supports all workflows (local, cloud, hybrid, evaluate) and execution profiles
# (wits, standard, docker, singularity, apptainer).
#
# Modes:
#   1. Interactive Guided Wizard (default when run without arguments)
#   2. Headless CLI Flag Mode (ideal for automated batch scripts and CI)
#
# Key Features:
#   - Universal Workflow Support: local, cloud (ICAv2 DRAGEN), hybrid, evaluate
#   - Profile-Agnostic: works across cluster (wits, slurm) and local/container engines
#   - Reference Harmonization: handles hg38 and GRCh38 builds seamlessly
#   - Auto-Generated Samplesheets & Biological Sex Info (via chrX/chrY coverage)
#   - Pre-flight Diagnostic Auditing before job submission
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Color formatting
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

# Defaults
INPUT_DIR=""
OUTDIR=""
WORKFLOW="local"
PROFILE="balanced"
EXECUTOR_PROFILE=""
SEXINFO_FILE=""
REF_DIR=""
GENOME_BUILD="GRCh38"
GENERATE_IGV=false
EXEC_MODE=""

# Auto-detect executor environment default
if [ -d "/dataG/ddd" ] && command -v sbatch &>/dev/null; then
    DEFAULT_EXECUTOR="wits"
elif command -v apptainer &>/dev/null; then
    DEFAULT_EXECUTOR="apptainer"
elif command -v singularity &>/dev/null; then
    DEFAULT_EXECUTOR="singularity"
elif command -v docker &>/dev/null; then
    DEFAULT_EXECUTOR="docker"
else
    DEFAULT_EXECUTOR="standard"
fi

print_help() {
    cat << HELP_EOF
CNV-Caller-Combiner Unified Launcher
Usage: ./launch.sh [OPTIONS]

Options:
  -i, --input-dir DIR          Directory containing cohort BAM or CRAM alignment files
  -o, --output-dir DIR         Destination directory for results and intermediate configs
  -w, --workflow MODE          Workflow mode: local (default), cloud, hybrid, evaluate
  -p, --profile PROFILE        Quality profile: balanced (default), high_precision, max_sensitivity
  -e, --executor-profile PROF  Environment profile: wits (default on cluster), standard, docker, singularity, apptainer
  -b, --build BUILD            Genome build: GRCh38 (default) or hg38
  -s, --sexinfo FILE           Custom sample sex metadata TSV (sample_id, sex). If omitted, auto-inferred.
  -r, --ref-dir DIR            Base directory for references (default checks Wits cluster and ./references)
  --batch                      Submit pipeline directly to SLURM cluster scheduler
  --local                      Run Nextflow pipeline directly in foreground terminal
  --summary [DIR]              Display executive summary report for a completed run directory
  --results [DIR]              Alias for --summary
  --enable-igv                 Generate standalone offline IGV HTML reports per sample
  --validate-only              Perform pre-flight diagnostics and generate configs without running
  -h, --help                   Show this help message and exit

Examples:
  Interactive wizard:
    ./launch.sh

  Headless SLURM submission (balanced local workflow):
    ./launch.sh -i /spaces/${USER}/cohort_bams -o /spaces/${USER}/cnv_run -p balanced --batch

  Inspect results from a completed run:
    ./launch.sh --summary /spaces/${USER}/cnv_run/balanced
HELP_EOF
    exit 0
}

# Parse CLI arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir) INPUT_DIR="$2"; shift 2 ;;
        -o|--output-dir) OUTDIR="$2"; shift 2 ;;
        -w|--workflow) WORKFLOW="$2"; shift 2 ;;
        -p|--profile|--filter-profile) PROFILE="$2"; shift 2 ;;
        -e|--executor-profile) EXECUTOR_PROFILE="$2"; shift 2 ;;
        -b|--build) GENOME_BUILD="$2"; shift 2 ;;
        -s|--sexinfo) SEXINFO_FILE="$2"; shift 2 ;;
        -r|--ref-dir) REF_DIR="$2"; shift 2 ;;
        --summary|--results)
            RESULTS_TARGET="${2:-}"
            if [ -z "${RESULTS_TARGET}" ] || [[ "${RESULTS_TARGET}" == -* ]]; then
                read -r -p "Enter path to results directory: " RESULTS_TARGET
            fi
            python3 "${SCRIPT_DIR}/bin/summarize_results.py" --results-dir "${RESULTS_TARGET}"
            exit 0
            ;;
        --batch) EXEC_MODE="batch"; shift ;;
        --local) EXEC_MODE="local"; shift ;;
        --validate-only) EXEC_MODE="validate_only"; shift ;;
        --enable-igv) GENERATE_IGV=true; shift ;;
        --skip-igv) GENERATE_IGV=false; shift ;;
        -h|--help) print_help ;;
        *) echo -e "${RED}Unknown option: $1${NC}" >&2; exit 1 ;;
    esac
done

EXECUTOR_PROFILE="${EXECUTOR_PROFILE:-${DEFAULT_EXECUTOR}}"

echo -e "${CYAN}${BOLD}"
echo "========================================================================"
echo "    🧬  CNV-CALLER-COMBINER: UNIFIED PIPELINE LAUNCHER"
echo "========================================================================"
echo -e "${NC}"

# 1. Determine workflow mode
if [ -z "${EXEC_MODE}" ] && [ -z "${INPUT_DIR}" ]; then
    echo -e "${BOLD}Step 1: Pipeline Workflow Selection${NC}"
    echo "1) local   [Default] - Run 6 local callers (GATK gCNV, XHMM, CANOES, CLAMMS, CNVkit, ExomeDepth, SvABA)"
    echo "2) cloud             - Run cloud-based Illumina ICAv2 DRAGEN Germline CNV/SV pipeline"
    echo "3) hybrid            - Run both local callers + ICAv2 DRAGEN cloud and merge all callsets"
    echo "4) evaluate          - Benchmark existing callsets against a truth BED"
    echo "5) summary           - View and summarize results from a completed run"
    read -r -p "Select workflow [1-5, default 1]: " WF_CHOICE
    case "${WF_CHOICE:-1}" in
        2) WORKFLOW="cloud" ;;
        3) WORKFLOW="hybrid" ;;
        4) WORKFLOW="evaluate" ;;
        5)
            read -r -p "Enter path to results directory: " RESULTS_TARGET
            python3 "${SCRIPT_DIR}/bin/summarize_results.py" --results-dir "${RESULTS_TARGET}"
            exit 0
            ;;
        *) WORKFLOW="local" ;;
    esac
fi
echo -e "Selected workflow: ${GREEN}${BOLD}${WORKFLOW}${NC}"

# 2. Input directory
if [ -z "${INPUT_DIR}" ]; then
    echo
    echo -e "${BOLD}Step 2: Input Sequencing Alignments${NC}"
    echo "Enter folder containing cohort BAM or CRAM files:"
    read -r -p "Folder path: " INPUT_DIR
fi

if [ ! -d "${INPUT_DIR}" ]; then
    echo -e "${RED}❌ Error: Directory '${INPUT_DIR}' does not exist.${NC}" >&2
    exit 1
fi
INPUT_DIR="$(cd "${INPUT_DIR}" && pwd)"

# 3. Output directory
if [ -z "${OUTDIR}" ]; then
    echo
    echo -e "${BOLD}Step 3: Results Destination${NC}"
    DEFAULT_OUT="/spaces/${USER}/cnv_${WORKFLOW}_$(date +%Y%m%d_%H%M)"
    read -r -p "Output folder [Default: ${DEFAULT_OUT}]: " OUTDIR
    OUTDIR="${OUTDIR:-${DEFAULT_OUT}}"
fi
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

# 4. Create samplesheet
echo
echo -e "${BOLD}Step 4: Creating and Validating Samplesheet...${NC}"
SAMPLESHEET="${OUTDIR}/samplesheet.tsv"
python3 "${SCRIPT_DIR}/bin/make_samplesheet.py" \
    --input-dir "${INPUT_DIR}" \
    --output "${SAMPLESHEET}"

TOTAL_SAMPLES=$(tail -n +2 "${SAMPLESHEET}" | wc -l)
if [ "${TOTAL_SAMPLES}" -eq 0 ]; then
    echo -e "${RED}❌ Error: No valid BAM/CRAM alignment pairs found in ${INPUT_DIR}.${NC}" >&2
    exit 1
fi
echo -e "Found ${GREEN}${TOTAL_SAMPLES} alignment samples${NC} in ${INPUT_DIR}"

IS_CRAM=false
if head -n 1 "${SAMPLESHEET}" | grep -qi "CRAM"; then
    IS_CRAM=true
fi

# 5. Sexinfo handling
echo
echo -e "${BOLD}Step 5: Sample Biological Sex Metadata${NC}"
if [ -n "${SEXINFO_FILE}" ] && [ -f "${SEXINFO_FILE}" ]; then
    FINAL_SEXINFO="${SEXINFO_FILE}"
    echo -e "Using provided sexinfo file: ${GREEN}${FINAL_SEXINFO}${NC}"
else
    FINAL_SEXINFO="${OUTDIR}/sexinfo.tsv"
    echo "Auto-inferring sample biological sex from chrX/chrY coverage via bin/make_sexinfo.py..."
    python3 "${SCRIPT_DIR}/bin/make_sexinfo.py" \
        --samplesheet "${SAMPLESHEET}" \
        --output "${FINAL_SEXINFO}" \
        --default-sex "F"
    echo -e "Generated sex metadata: ${GREEN}${FINAL_SEXINFO}${NC}"
fi

# 6. Quality Profile Selection
if [ -z "${EXEC_MODE}" ]; then
    echo
    echo -e "${BOLD}Step 6: Quality Filtering Profile${NC}"
    echo "1) Balanced [Recommended]   - Optimum sensitivity and precision (F1-score maximized)"
    echo "2) High Precision           - Minimal false positives (strict caller intersection)"
    echo "3) Max Sensitivity          - Maximum discovery of subtle/low-frequency CNVs"
    read -r -p "Select [1-3, default 1]: " PROFILE_CHOICE
    PROFILE_CHOICE="${PROFILE_CHOICE:-1}"
    case "${PROFILE_CHOICE}" in
        2) PROFILE="high_precision" ;;
        3) PROFILE="max_sensitivity" ;;
        *) PROFILE="balanced" ;;
    esac
fi
echo -e "Active profile: ${GREEN}${BOLD}${PROFILE}${NC}"

# 7. Reference Bundle Discovery (Supporting both hg38 and GRCh38)
echo
echo -e "${BOLD}Step 7: Reference Discovery (${GENOME_BUILD})...${NC}"
WITS_REF="/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
LOCAL_REF="${SCRIPT_DIR}/references/fasta/Homo_sapiens_assembly38.fasta"

if [ -n "${REF_DIR}" ] && [ -d "${REF_DIR}" ]; then
    REF_FASTA=$(find "${REF_DIR}" -maxdepth 2 -name "*.fa" -o -name "*.fasta" | head -n 1)
elif [ -f "${WITS_REF}" ]; then
    REF_FASTA="${WITS_REF}"
elif [ -f "${LOCAL_REF}" ]; then
    REF_FASTA="${LOCAL_REF}"
else
    REF_FASTA=""
fi

if [ -z "${REF_FASTA}" ] || [ ! -f "${REF_FASTA}" ]; then
    echo -e "${YELLOW}⚠️  Standard reference FASTA not detected automatically.${NC}"
    echo "Run '${BOLD}./bin/setup_reference_bundle.sh${NC}' to stage standard reference assets."
    if [ -z "${EXEC_MODE}" ]; then
        read -r -p "Enter path to reference FASTA: " REF_FASTA
    else
        echo -e "${RED}❌ Missing required reference FASTA. Provide via --ref-dir or stage references first.${NC}" >&2
        exit 1
    fi
fi

REF_FAI="${REF_FASTA}.fai"
REF_DICT="${REF_FASTA%.*}.dict"

TARGETS_BED="/dataG/ddd/data/resources/canoes/probes_sanger.bed"
if [ ! -f "${TARGETS_BED}" ]; then
    TARGETS_BED="${SCRIPT_DIR}/references/targets.bed"
fi

CLAMMS_MAP="/dataG/ddd/data/resources/clamms/mappability.bed"
CLAMMS_SPEC="/dataG/ddd/data/resources/clamms/clamms_special_regions.grch38_sorted.bed"
if [ ! -f "${CLAMMS_MAP}" ]; then
    CLAMMS_MAP="${SCRIPT_DIR}/references/clamms/mappability.bed"
    CLAMMS_SPEC="${SCRIPT_DIR}/references/clamms/clamms_special_regions.bed"
fi

XHMM_CONF="/dataG/ddd/data/resources/xhmm/params.txt"
if [ ! -f "${XHMM_CONF}" ]; then
    XHMM_CONF="${SCRIPT_DIR}/references/xhmm/params.txt"
fi

ANNOTSV_DIR="${SCRIPT_DIR}/annotsv_db/AnnotSV_annotations/share/AnnotSV"

# 8. Generate params JSON
PARAMS_FILE="${OUTDIR}/params-${WORKFLOW}-${PROFILE}.json"
echo
echo -e "${BOLD}Step 8: Generating Configuration File...${NC}"

cat << JSON_EOF > "${PARAMS_FILE}"
{
    "workflow": "${WORKFLOW}",
    "outdir": "${OUTDIR}/${PROFILE}",
    "bind_paths": "/spaces,/dataG,/home/${USER},/dataB/aux,/dataG/ddd,/dataG/ddd-2023",
    "genome_build": "${GENOME_BUILD}",
    "build": "${GENOME_BUILD}",
    "runOptions": "",
    
    $( [ "${IS_CRAM}" = true ] && echo "\"samplesheet_crams\": \"${SAMPLESHEET}\"," || echo "\"samplesheet_bams\": \"${SAMPLESHEET}\"," )
    "ref": "${REF_FASTA}",
    "fai": "${REF_FAI}",
    "fasta": "${REF_FASTA}",
    "dict": "${REF_DICT}",
    
    "targets": "${TARGETS_BED}",
    "probes": "${TARGETS_BED}",
    
    "mappability": "${CLAMMS_MAP}",
    "special_reg": "${CLAMMS_SPEC}",
    "sexinfo": "${FINAL_SEXINFO}",
    
    "xhmm_conf": "${XHMM_CONF}",
    "annotsv_db": "${ANNOTSV_DIR}",
    "annotsv_version": "3.5.10",
    
    "filter_profile": "${PROFILE}",
    "merger_mode": "truvari",
    "run_svaba": true,
    "run_refinement": true,
    "generate_igv_reports": ${GENERATE_IGV},
    "skip_cohort_merge": false,
    "skip_af_filtering": false,
    "skip_annotation": false,
    "skip_ml_training": true
}
JSON_EOF
echo -e "Saved parameters: ${GREEN}${PARAMS_FILE}${NC}"

# 9. Pre-flight diagnostics
echo
echo -e "${BOLD}Step 9: Running Pre-flight Diagnostic Audit...${NC}"
python3 "${SCRIPT_DIR}/bin/validate_run.py" --params "${PARAMS_FILE}" || {
    echo -e "${YELLOW}⚠️  Pre-flight validation reported warnings. Review above messages.${NC}"
}

# 10. Execution dispatch
if [ -z "${EXEC_MODE}" ]; then
    echo
    echo -e "${BOLD}Step 10: Execution Mode${NC}"
    echo "1) Submit to SLURM Cluster (Recommended for cohorts) [Default]"
    echo "2) Run directly in terminal (Foreground, for small test batches)"
    echo "3) Validate only (Save configs and exit)"
    read -r -p "Select [1-3, default 1]: " LAUNCH_CHOICE
    case "${LAUNCH_CHOICE:-1}" in
        2) EXEC_MODE="local" ;;
        3) EXEC_MODE="validate_only" ;;
        *) EXEC_MODE="batch" ;;
    esac
fi

load_environment() {
    if command -v module &>/dev/null; then
        module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7 2>/dev/null || true
    fi
    if ! command -v nextflow &>/dev/null; then
        echo -e "${RED}❌ Nextflow is not installed or not in PATH.${NC}" >&2
        echo "Install Nextflow (curl -s https://get.nextflow.io | bash) or load its environment module." >&2
        exit 1
    fi
}

case "${EXEC_MODE}" in
    validate_only)
        echo
        echo -e "${GREEN}Configuration and validation complete!${NC}"
        echo -e "  Parameters:  ${PARAMS_FILE}"
        echo -e "  Samplesheet: ${SAMPLESHEET}"
        echo -e "  Sex Info:    ${FINAL_SEXINFO}"
        echo "To execute manually:"
        echo "  nextflow run cnv-caller-combiner.nf -profile ${EXECUTOR_PROFILE},${PROFILE} -params-file ${PARAMS_FILE} -resume"
        exit 0
        ;;

    local)
        echo
        echo -e "${CYAN}🚀 Launching Nextflow pipeline in foreground...${NC}"
        load_environment
        nextflow run "${SCRIPT_DIR}/cnv-caller-combiner.nf" \
            -profile "${EXECUTOR_PROFILE},${PROFILE}" \
            -params-file "${PARAMS_FILE}" \
            --workflow "${WORKFLOW}" \
            --outdir "${OUTDIR}/${PROFILE}" \
            -w "${OUTDIR}/work" \
            -resume
        ;;

    batch|*)
        SBATCH_SCRIPT="${OUTDIR}/run_pipeline.sbatch"
        cat << SBATCH_EOF > "${SBATCH_SCRIPT}"
#!/usr/bin/env bash
#SBATCH --job-name=cnv-${WORKFLOW}
#SBATCH --output=${OUTDIR}/slurm-%x-%j.out
#SBATCH --error=${OUTDIR}/slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail
if command -v module &>/dev/null; then
    module load java/21 nextflow/25.10 samtools/1.23 bedtools/2.30.0 SURVIVOR/1.0.7 2>/dev/null || true
fi

export NXF_SINGULARITY_CACHEDIR="\${NXF_SINGULARITY_CACHEDIR:-/spaces/\${USER}/nxf_cache}"
export NXF_ANSI_LOG=false

cd "${SCRIPT_DIR}"

nextflow run cnv-caller-combiner.nf \\
    -profile ${EXECUTOR_PROFILE},medium,${PROFILE} \\
    -params-file "${PARAMS_FILE}" \\
    --workflow "${WORKFLOW}" \\
    --outdir "${OUTDIR}/${PROFILE}" \\
    -w "${OUTDIR}/work" \\
    -resume

echo "Pipeline completed successfully at \$(date)"
if [ -f "${SCRIPT_DIR}/bin/summarize_results.py" ]; then
    python3 "${SCRIPT_DIR}/bin/summarize_results.py" --results-dir "${OUTDIR}/${PROFILE}" || true
fi
SBATCH_EOF
        chmod +x "${SBATCH_SCRIPT}"

        echo
        echo -e "${CYAN}Submitting job to SLURM cluster scheduler...${NC}"
        JOB_SUB=$(sbatch "${SBATCH_SCRIPT}")
        JOB_ID=$(echo "${JOB_SUB}" | awk '{print $NF}')
        echo -e "${GREEN}${BOLD}🎉 Job successfully submitted! (Job ID: ${JOB_ID})${NC}"
        echo
        echo "Helpful monitoring commands:"
        echo -e "  Check job queue:   ${BOLD}squeue -j ${JOB_ID}${NC}"
        echo -e "  Follow live logs:  ${BOLD}tail -f ${OUTDIR}/slurm-cnv-${WORKFLOW}-${JOB_ID}.out${NC}"
        echo -e "  Summarize results: ${BOLD}python3 ${SCRIPT_DIR}/bin/summarize_results.py --results-dir ${OUTDIR}/${PROFILE}${NC}"
        echo
        ;;
esac
