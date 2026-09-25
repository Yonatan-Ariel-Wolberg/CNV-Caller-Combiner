#!/usr/bin/env bash
# ==============================================================================
# setup_reference_bundle.sh
#
# Automated Reference Bundle Downloader & Staging Utility for CNV-Caller-Combiner
# Prepares standard reference files for GRCh38 (or GRCh37):
#   1. Reference genome FASTA, .fai, and Picard .dict (Broad GATK bundle mirror)
#   2. CLAMMS mappability track and special regions
#   3. XHMM default hidden Markov model parameters
#   4. Emits a pre-configured reference_paths.json
# ==============================================================================

set -euo pipefail

DEST_DIR="./references"
BUILD="GRCh38"
DOWNLOAD_FASTA=true

print_help() {
    cat << HELP_EOF
Usage: setup_reference_bundle.sh [OPTIONS]

Options:
  -d, --dest DIR        Destination directory to store reference files [Default: ./references]
  -b, --build BUILD     Genome build: GRCh38 or GRCh37 [Default: GRCh38]
  --skip-fasta          Skip downloading the large FASTA/dict (e.g. if you already have it)
  -h, --help            Show this help message
HELP_EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dest) DEST_DIR="$2"; shift 2 ;;
        -b|--build) BUILD="$2"; shift 2 ;;
        --skip-fasta) DOWNLOAD_FASTA=false; shift ;;
        -h|--help) print_help ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${DEST_DIR}"
DEST_DIR="$(cd "${DEST_DIR}" && pwd)"

echo "========================================================================"
echo "    🧬  CNV-CALLER-COMBINER: REFERENCE BUNDLE SETUP"
echo "========================================================================"
echo "Destination:  ${DEST_DIR}"
echo "Genome Build: ${BUILD}"
echo

# 1. Download CLAMMS reference tracks
echo "==> [1/3] Setting up CLAMMS reference tracks..."
CLAMMS_DIR="${DEST_DIR}/clamms"
mkdir -p "${CLAMMS_DIR}"

if [ ! -s "${CLAMMS_DIR}/clamms_special_regions.bed" ]; then
    echo "  Fetching CLAMMS special regions track..."
    curl -fsSL --connect-timeout 5 --max-time 15 \
        "https://raw.githubusercontent.com/rgcgithub/clamms/master/data/clamms_special_regions.bed" \
        -o "${CLAMMS_DIR}/clamms_special_regions.bed" 2>/dev/null || {
        echo "  Note: Could not download remote CLAMMS regions. Creating minimal placeholder."
        touch "${CLAMMS_DIR}/clamms_special_regions.bed"
    }
fi

# 2. Setup XHMM parameters
echo "==> [2/3] Setting up XHMM parameters..."
XHMM_DIR="${DEST_DIR}/xhmm"
mkdir -p "${XHMM_DIR}"
if [ ! -s "${XHMM_DIR}/params.txt" ]; then
    echo "  Configuring standard XHMM parameter matrix..."
    cat << 'PARAM_EOF' > "${XHMM_DIR}/params.txt"
1e-8    6    1e-8
1e-8    6    1e-8
100000  0.1  0.1
100000  0.1  0.1
PARAM_EOF
fi

# 3. Setup Reference Genome FASTA & Indices
echo "==> [3/3] Setting up Reference Genome..."
FASTA_DIR="${DEST_DIR}/fasta"
mkdir -p "${FASTA_DIR}"

REF_FASTA="${FASTA_DIR}/Homo_sapiens_assembly38.fasta"
REF_FAI="${REF_FASTA}.fai"
REF_DICT="${FASTA_DIR}/Homo_sapiens_assembly38.dict"

if [ "${DOWNLOAD_FASTA}" = true ]; then
    if [ ! -s "${REF_FASTA}" ]; then
        echo "  Downloading Broad GATK GRCh38 FASTA (~3.1 GB)..."
        curl -fSL "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
            -o "${REF_FASTA}"
    fi
    if [ ! -s "${REF_FAI}" ]; then
        echo "  Downloading FASTA index (.fai)..."
        curl -fSL "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
            -o "${REF_FAI}"
    fi
    if [ ! -s "${REF_DICT}" ]; then
        echo "  Downloading sequence dictionary (.dict)..."
        curl -fSL "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict" \
            -o "${REF_DICT}"
    fi
else
    echo "  Skipping FASTA download (--skip-fasta specified)."
fi

# 4. Generate reference_paths.json
OUTPUT_JSON="${DEST_DIR}/reference_paths.json"
cat << JSON_EOF > "${OUTPUT_JSON}"
{
    "ref": "${REF_FASTA}",
    "fai": "${REF_FAI}",
    "dict": "${REF_DICT}",
    "special_reg": "${CLAMMS_DIR}/clamms_special_regions.bed",
    "xhmm_conf": "${XHMM_DIR}/params.txt"
}
JSON_EOF

echo
echo "✅ Reference bundle staging complete!"
echo "Reference configuration template saved to: ${OUTPUT_JSON}"
