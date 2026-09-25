# XHMM (eXome Hidden Markov Model)

## 1. Methodology & Algorithm
**XHMM** is a gold-standard statistical tool designed to discover and genotype germline copy-number variation from whole-exome sequencing read depth. Its core breakthrough is using **Principal Component Analysis (PCA)** to remove systematic technical noise across capture batches:
* **Target Interval Depth Matrix**: Constructs a unified $N \times M$ matrix representing the depth of $N$ samples across $M$ targeted intervals.
* **PCA Noise Removal**: Deconstructs the centered and standardized read-depth matrix using Singular Value Decomposition (SVD). High-variance components that reflect technical noise (GC content, bait efficiency, hybridization temperature gradients) are projected out.
* **Matrix Normalization & Z-Scoring**: Normalized read depths are converted into interval-specific z-scores centered around cohort expectations.
* **Viterbi HMM Decoding & Genotyping**: Employs an HMM with three underlying copy number states:
  - Deletion (copy number 1 or 0)
  - Diploid normal (copy number 2)
  - Duplication (copy number 3 or higher)
  Genotyping calculates the Phred-scaled likelihood (`SQ` score) that a sample harbors a non-diploid event over each segment.

---

## 2. Nextflow Implementation & Processes
The XHMM module is encapsulated in `modules/callers/modules-xhmm.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `GATK_DEPTH_OF_COVERAGE` | `gatk` | Computes per-target read counts across all cohort BAMs using GATK `DepthOfCoverage`. |
| `MERGE_DEPTH_MATRICES` | `xhmm` | Compiles individual sample depth tables into a single matrix across all targeted intervals. |
| `XHMM_FILTER_TARGETS` | `xhmm` | Removes intervals with extreme GC content, low mean depth (< 10x), or excessive variance. |
| `XHMM_PCA_NORMALIZATION` | `xhmm` | Performs SVD and projects out principal components explaining technical cohort-wide variance. |
| `XHMM_HMM_DISCOVERY` | `xhmm` | Runs Viterbi decoding over normalized z-score matrices to discover CNV intervals. |
| `XHMM_GENOTYPING` | `xhmm` | Generates posterior genotype probabilities and Phred-scaled quality scores (`SQ`, `NQ`). |
| `FORMAT_XHMM_VCF` | `python` | Converts native XHMM VCF output to standardized VCF 4.2 with canonical `SVLEN` and `END` tags. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA file and `.fai` index.
* **Configuration Parameters**: `params/xhmm.conf` controlling PCA variance filtering thresholds.

---

## 4. Key Parameters & Quality Thresholds
* `params.filter_xhmm_qscore`: Quality score threshold for `SQ` (default: `>= 60` in balanced mode).
* `params.filter_xhmm_probes`: Minimum number of supporting intervals (default: `>= 3`).
* `params.filter_xhmm_min_size`: Minimum length in base pairs (default: `1000 bp`).

---

## 5. Output Files
* `out_XHMM/raw/<SAMPLE_ID>.xhmm.raw.vcf.gz`: Raw output VCF containing candidate calls.
* `out_XHMM/filtered/<SAMPLE_ID>.xhmm.filtered.vcf.gz`: Quality-filtered and normalized callset.
