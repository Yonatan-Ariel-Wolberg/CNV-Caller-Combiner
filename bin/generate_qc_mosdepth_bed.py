#!/usr/bin/env python3
import os
import glob
import subprocess

def main():
    failed_file = "/home/ywolberg/CNV-Caller-Combiner/qc_results_ddd_uk/summary_report/failed_samples.txt"
    with open(failed_file) as f:
        failed_samples = set(line.strip().replace('.md', '') for line in f if line.strip())

    print(f"Loaded {len(failed_samples)} failed samples.")

    consensus_dir = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/100x-mosdepth/consensus"
    all_bed_files = glob.glob(os.path.join(consensus_dir, "*/*.usable_20x.bed"))
    print(f"Found {len(all_bed_files)} total usable_20x.bed files.")

    qc_passed_beds = []
    for bed_path in all_bed_files:
        filename = os.path.basename(bed_path)
        sample_id = filename.split('.')[0]
        if sample_id not in failed_samples:
            qc_passed_beds.append(bed_path)

    print(f"Collected {len(qc_passed_beds)} QC-passed BED files.")
    assert len(qc_passed_beds) == 63, f"Expected 63 QC-passed BED files, but found {len(qc_passed_beds)}"

    output_dir = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/mosdepth"
    os.makedirs(output_dir, exist_ok=True)
    temp_intersection = os.path.join(output_dir, "qc_passed_intersection.bed")
    final_bed = os.path.join(output_dir, "mosdepth_sanger_qc_passed.bed")

    # Run bedtools multiinter via Apptainer
    print("Running bedtools multiinter inside Apptainer...")
    container = "/home/ywolberg/.apptainer/cache/quay.io-biocontainers-bedtools-2.31.1--h13024bc_3.img"
    cmd = [
        "apptainer", "exec",
        "-B", "/spaces,/dataG,/home/ywolberg",
        container,
        "bedtools", "multiinter", "-i"
    ] + qc_passed_beds
    with open(temp_intersection, "w") as out_fh:
        subprocess.run(cmd, stdout=out_fh, check=True)

    # Filter for regions where >= 90% (57 of 63) samples passed
    print("Filtering consensus regions (>=57 samples)...")
    count_kept = 0
    with open(temp_intersection) as in_fh, open(final_bed, "w") as out_fh:
        for line in in_fh:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                nsamples = int(parts[3])
                if nsamples >= 57:
                    out_fh.write(line)
                    count_kept += 1

    print(f"Consensus core BED created at {final_bed} with {count_kept} regions.")
    # Clean up temp file
    if os.path.exists(temp_intersection):
        os.remove(temp_intersection)

if __name__ == "__main__":
    main()
