#!/bin/bash
#SBATCH --job-name=cnv_combiner_multi
#SBATCH --output=run_local_combinations.out
#SBATCH --error=run_local_combinations.err
#SBATCH --time=96:00:00
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Define paths
CONTROL_BED="/dataG/ddd/data/resources/canoes/probes_sanger.bed"
MOSDEPTH_BED="/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/mosdepth/mosdepth_sanger.bed"
MOSDEPTH_BED_QC_PASSED="/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/mosdepth/mosdepth_sanger_qc_passed.bed"
PARAMS_FILE="/home/ywolberg/CNV-Caller-Combiner/params/ddd-uk/params-local-master.json"

# Configurations to run:
# 1. QC=false, BED=Control, Merger=SURVIVOR
echo "Running combination 1: QC=false, BED=Control, Merger=SURVIVOR"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc false \
  --probes_bed $CONTROL_BED \
  --merger_mode survivor \
  -resume

# 2. QC=false, BED=Control, Merger=Truvari
echo "Running combination 2: QC=false, BED=Control, Merger=Truvari"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc false \
  --probes_bed $CONTROL_BED \
  --merger_mode truvari \
  -resume

# 3. QC=true, BED=Control, Merger=SURVIVOR
echo "Running combination 3: QC=true, BED=Control, Merger=SURVIVOR"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc true \
  --probes_bed $CONTROL_BED \
  --merger_mode survivor \
  -resume

# 4. QC=true, BED=Control, Merger=Truvari
echo "Running combination 4: QC=true, BED=Control, Merger=Truvari"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc true \
  --probes_bed $CONTROL_BED \
  --merger_mode truvari \
  -resume

# 5. QC=false, BED=Mosdepth, Merger=SURVIVOR
echo "Running combination 5: QC=false, BED=Mosdepth, Merger=SURVIVOR"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc false \
  --probes_bed $MOSDEPTH_BED \
  --merger_mode survivor \
  -resume

# 6. QC=false, BED=Mosdepth, Merger=Truvari
echo "Running combination 6: QC=false, BED=Mosdepth, Merger=Truvari"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc false \
  --probes_bed $MOSDEPTH_BED \
  --merger_mode truvari \
  -resume

# 7. QC=true, BED=Mosdepth, Merger=SURVIVOR
echo "Running combination 7: QC=true, BED=Mosdepth, Merger=SURVIVOR"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc true \
  --probes_bed $MOSDEPTH_BED_QC_PASSED \
  --merger_mode survivor \
  -resume

# 8. QC=true, BED=Mosdepth, Merger=Truvari
echo "Running combination 8: QC=true, BED=Mosdepth, Merger=Truvari"
nextflow run main.nf -profile wits,small -params-file $PARAMS_FILE \
  --workflow local \
  --run_qc true \
  --probes_bed $MOSDEPTH_BED_QC_PASSED \
  --merger_mode truvari \
  -resume
