# 🧬 Beginner's Guide to CNV-Caller-Combiner

Welcome! If you are a researcher, clinician, or student with limited bioinformatics or programming experience, this guide is written for you. 

The **CNV-Caller-Combiner** automatically discovers **Copy Number Variations (CNVs)**—deletions and duplications of genomic segments—from Whole Exome Sequencing (WES) data.

---

## 💡 Why Use This Pipeline?

Single CNV callers frequently suffer from high false-positive rates or miss real variants due to PCR bias, GC fluctuations, and uneven capture efficiency:
* **CANOES** & **CLAMMS** excel at read-depth modeling across targeted intervals.
* **XHMM** uses Principal Component Analysis (PCA) to remove systematic cohort noise.
* **CNVkit** models both on-target and off-target read depths.
* **GATK gCNV** applies Bayesian probabilistic models across intervals.
* **ExomeDepth** dynamically selects customized reference cohorts per sample.
* **SvABA** reconstructs structural variant (SV) breakpoints at single-base-pair resolution.

**This pipeline runs these tools together, merges their consensus using geometric overlap, refines the breakpoints to nucleotide resolution using size-aware constraints, annotates with AnnotSV, and classifies variants into clinical ACMG Pathogenic tiers—all automatically.**

---

## 🚀 The Quickstart Launcher (`launch.sh`)

You do not need to write complex command-line scripts or edit configuration files by hand.

### Mode A: Interactive Guided Wizard (Zero Coding)

Navigate to the repository folder and run:
```bash
cd /home/ywolberg/CNV-Caller-Combiner
./launch.sh
```

The launcher will walk you through a setup wizard:
1. **Workflow Selection**: Choose `local` (default: 6 callers), `cloud` (ICAv2 DRAGEN), `hybrid` (both), or `evaluate`.
2. **Input alignments folder**: Paste the path to your folder containing `.bam` or `.cram` files (e.g. `/dataG/my_cohort/`).
3. **Results folder**: Enter where you want your output saved (or press `Enter` to use the default).
4. **Sex Information**: Automatically inferred from chrX / chrY coverage using `bin/make_sexinfo.py`, or provide a custom file.
5. **Filtering Profile**:
   * Press `1` for **Balanced** *(Recommended for most studies; F1-score maximized)*.
   * Press `2` for **High Precision** *(Strict thresholds, minimal false positives; ideal if validating with qPCR)*.
   * Press `3` for **Max Sensitivity** *(Relaxed thresholds; ideal for undiagnosed rare diseases)*.
6. **Launch Mode**: Press `1` to submit the job to the cluster scheduler (SLURM) in the background.

### Mode B: Headless Command-Line Flags (For Advanced Users & Batch Scripts)

You can bypass prompts completely by passing flags:
```bash
# Submit local workflow to SLURM cluster:
./launch.sh -i /spaces/$USER/cohort_bams -o /spaces/$USER/cnv_run -p balanced --batch

# Run directly in foreground terminal:
./launch.sh -i /spaces/$USER/cohort_bams -o /spaces/$USER/cnv_run -p balanced --local

# Pre-flight check only (verify files without running):
./launch.sh -i /spaces/$USER/cohort_bams -o /spaces/$USER/cnv_run --validate-only
```

---

## 🔍 How to Monitor Your Run

You can check on your job anytime with these commands:

```bash
# Check if your job is running or queued:
squeue -u $USER

# Watch live pipeline output:
tail -f /path/to/results/slurm-cnv-local-<JOB_ID>.out
```

---

## 📊 How to Read Your Results (Without Coding!)

When the pipeline completes, navigate to your results folder. You will find:

### 1. The Executive Clinical Summary (`CLINICAL_SUMMARY_REPORT.md`)
Run the summary tool anytime to get an instant breakdown of findings:
```bash
python3 bin/summarize_results.py --results-dir /path/to/results/balanced
```
It displays:
* Total Pathogenic (`ACMG_5`) and Likely Pathogenic (`ACMG_4`) calls.
* Top recurrent disease genes affected across your cohort.
* Direct paths to tabular files ready for spreadsheet analysis.

### 2. Opening in Microsoft Excel or Google Sheets
Open the file:
```
consensus_evaluate/*ACMG_4_5.tsv
```
This spreadsheet-ready table contains high-confidence clinical candidates with the following columns:
* **`SampleID`**: The patient / sample identifier.
* **`CHROM`**, **`POS`**, **`END`**: Refined genomic coordinates (GRCh38).
* **`SVTYPE`**: `DEL` (Deletion) or `DUP` (Duplication).
* **`CALLERS`**: Which callers agreed on this variant (e.g. `CANOES;CLAMMS;XHMM;GCNV`).
* **`NUM_CALLERS`**: Total number of callers supporting the call.
* **`ACMG_CLASS`**: `5` (Pathogenic) or `4` (Likely Pathogenic).
* **`GENES`**: Specific disease-associated genes overlapped by the CNV.
* **`MUTATIONAL_MECHANISM`**: Biophysical formation mechanism (e.g. `NHEJ`, `MMBIR_FOSTES`, `NAHR_CANDIDATE`).

### 3. Viewing in IGV (Integrative Genomics Viewer)
To visualize reads confirming a deletion or duplication:
1. Open IGV on your desktop.
2. Select genome **Human (hg38)**.
3. Load the patient BAM: `/path/to/results/converted_bams/<SampleID>.bam`.
4. Load the consensus VCF: `/path/to/results/consensus_refined/<SampleID>.refined.vcf.gz`.
5. Jump to the coordinate (e.g. `chr7:117,550,000-117,600,000`) to visually confirm read-depth drop or duplication.

---

## 🛠️ Built-in Helper Tools for Beginners

| Tool | Command | What It Does |
| :--- | :--- | :--- |
| **Reference Staging Utility** | `./bin/setup_reference_bundle.sh -d ./references` | Downloads and stages CLAMMS, XHMM, and Broad GATK GRCh38 references in one command. |
| **Sex Metadata Inferrer** | `python3 bin/make_sexinfo.py -s samples.tsv -o sexinfo.tsv` | Automatically infers sample biological sex from chrX/chrY coverage. |
| **Samplesheet Maker** | `python3 bin/make_samplesheet.py -i /path/to/bams -o samples.tsv` | Automatically scans any folder, checks for `.bai`/`.crai` indexes, and generates a valid samplesheet. |
| **Pre-Flight Checker** | `python3 bin/validate_run.py -p params.json` | Runs a 10-point audit verifying all reference FASTA files, BEDs, and databases before you start. |
| **Results Summarizer** | `python3 bin/summarize_results.py -r /path/to/results` | Tallies ACMG classes and lists candidate genes in a clean summary table. |

---

## ❓ Frequently Asked Questions (FAQ)

#### Does this pipeline handle both `hg38` and `GRCh38`?
**Yes.** All analysis scripts use coordinate harmonization (`standardize_chrom`) to automatically recognize and harmonize both UCSC `chr1` and Ensembl `1` naming styles.

#### Does AnnotSV need to be installed manually?
**No.** If no local database is provided, the pipeline automatically downloads, verifies, compiles, and configures AnnotSV into `${projectDir}/annotsv_db` during the first run.

#### What if some of my BAMs are missing index (`.bai`) files?
Run `samtools index my_sample.bam` to generate the `.bai` index. The launcher and pre-flight audit will warn you if any indexes are missing or corrupted.

#### Can I re-run with a different filtering profile without re-running all the callers?
**Yes!** Nextflow's `-resume` feature caches caller results. If you change `-p balanced` to `-p high_precision`, it skips the time-consuming caller runs and re-filters in minutes.
