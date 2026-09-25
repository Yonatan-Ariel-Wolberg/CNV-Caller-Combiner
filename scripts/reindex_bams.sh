#!/usr/bin/env bash
# =============================================================================
# reindex_bams.sh: Re-indexes all BAMs in parallel to standard BAI format (BAI\x01)
# =============================================================================
set -euo pipefail

BAM_DIR="/spaces/ywolberg/2-test-cnv-caller-combiner/converted_bams"
module load samtools/1.23 || true

echo "Scanning BAM files in ${BAM_DIR}..."

reindex_one() {
    bam="$1"
    bai="${bam}.bai"
    magic=$(head -c 4 "$bai" 2>/dev/null || true)
    if [ "$magic" != $'BAI\x01' ]; then
        echo "[RE-INDEX] Re-indexing $(basename "$bam") to true BAI format..."
        rm -f "$bai"
        samtools index -b -@ 2 "$bam"
        touch -r "$bam" "$bai" 2>/dev/null || true
    else
        echo "[OK] $(basename "$bam") already has valid BAI format."
    fi
}
export -f reindex_one

find "${BAM_DIR}" -name "*.bam" | xargs -n 1 -P 12 -I {} bash -c 'reindex_one "{}"'

echo "=========================================================================="
echo "Verification Summary:"
python3 -c '
import glob
bais = glob.glob("/spaces/ywolberg/2-test-cnv-caller-combiner/converted_bams/*.bam.bai")
valid = 0
for f in bais:
    with open(f, "rb") as fh:
        if fh.read(4) == b"BAI\x01":
            valid += 1
print(f"Total BAI files: {len(bais)}")
print(f"Valid standard BAI (BAI\\x01): {valid}")
assert valid == len(bais) and len(bais) >= 90, "Not all BAMs have valid BAI indices!"
print("All BAM indices verified successfully!")
'
