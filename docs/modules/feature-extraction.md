# Feature Extraction Module

## 1. Methodology & Algorithm
The **Feature Extraction Module** translates raw genomic consensus calls and BAM alignments into rich numerical feature matrices for downstream Machine Learning classification (CatBoost / XGBoost):
* **Genomic Architecture Features**:
  - `SVLEN`: Absolute length of the structural event.
  - `GC_CONTENT`: Regional GC percentage across the variant span.
  - `NUM_TARGETS`: Number of targeted capture exons overlapped.
  - `TARGET_DENSITY`: Exons per kilobase.
* **Caller Concordance Features**:
  - `NUM_CALLERS`: Total count of individual callers supporting the call.
  - One-hot indicator variables for each caller (`HAS_CANOES`, `HAS_CLAMMS`, `HAS_XHMM`, `HAS_CNVKIT`, `HAS_GCNV`, `HAS_EXOMEDEPTH`, `HAS_DRAGEN`).
* **Depth & Quality Features**:
  - Normalized Phred-scaled caller quality scores (`QUAL_CANOES`, `QUAL_CLAMMS`, etc.).
  - Coverage ratio deviation: Observed vs expected read depth within the segment.
* **Cohort Context Features**:
  - `COHORT_AF`: Internal cohort allele frequency derived across the entire dataset.

---

## 2. Nextflow Implementation & Processes
The feature extraction module is encapsulated in `modules/ml/modules-feature-extraction.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `COMPUTE_TARGET_COVERAGE` | `mosdepth` | Rapidly computes per-base and per-target depth distributions across BAMs using `mosdepth`. |
| `FEATURE_EXTRACTION` | `python` | Extracts interval annotations, GC statistics, and caller metrics into per-sample tabular feature matrices (`.tsv`). |
| `BUILD_COHORT_FEATURE_DATABASE` | `python` | Merges per-sample TSVs and Mosdepth metrics into a unified cohort-wide feature matrix SQLite/Parquet database. |

---

## 3. Inputs and Reference Requirements
* **Consensus VCF**: Refined intra-sample consensus VCF.
* **Alignments / Depth Files**: Sample BAM/CRAM files or Mosdepth coverage summaries.
* **Reference Genome & BED**: Reference FASTA and target exon BED file.

---

## 4. Output Files
* `out_FEATURES/<SAMPLE_ID>_features.tsv`: Tabular feature matrix for the sample.
* `ml_features/cohort_features_database.parquet`: Unified feature matrix across all cohort samples.
