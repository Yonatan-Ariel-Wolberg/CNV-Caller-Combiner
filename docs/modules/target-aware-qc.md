# Target-Aware Quality Control (QC) & Sample Gating Module

## 1. Methodology & Algorithm
Exome sequencing datasets exhibit capture kit variability, uneven hybridization efficiency, and batch-specific library complexity. The **Target-Aware QC Module** provides robust quality assessment and sample gating prior to running heavy caller algorithms:
* **Target-Aware Coverage Breadth & Depth**: Rather than calculating broad whole-genome metrics, it restricts analysis strictly to the targeted capture intervals. It computes target breadth at $\ge 10\times$ and $\ge 20\times$ coverage.
* **Cohort Correlation Matrix**: Computes pairwise Pearson correlation coefficients of target depth across all samples. Samples exhibiting low correlation against the cohort median ($r < 0.85$ or $r < 0.90$) represent degraded DNA, capture failure, or severe batch anomalies.
* **Calibrated Sanger Exome Thresholds**: Standard exome QC standards (e.g. $\ge 90\%$ breadth at 10x) often spuriously flag older Sanger capture designs as failures due to narrower target enrichment profiles. This module provides specialized, validated thresholds for Sanger exomes ($75.0\%$ breadth at 10x, $0.85$ correlation).
* **FIPS Bypass Mount in Containers**: On modern enterprise Linux host nodes (such as Rocky Linux 9 or RHEL 9) where the Federal Information Processing Standard (FIPS) is enabled, OpenJDK inside unprivileged Apptainer containers often crashes due to crypto-module integrity mismatches. This module integrates an automatic bind-mount bypass (`/proc/sys/crypto/fips_enabled` zero-stub) to ensure container stability.
* **Dynamic 3-Column vs 4-Column BED Handling**: Automatically detects whether input target BED files contain 3 columns (chrom, start, end) or 4 columns (including target name), preventing coordinate-offset errors during `samtools bedcov` parsing.

---

## 2. Nextflow Implementation & Processes
Encapsulated in `modules/sample-qc/target-aware-exome-sequencing-quality-assess.nf`:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `SAMTOOLS_QC_BAM` / `SAMTOOLS_QC_CRAM` | `samtools` | Performs BAM integrity checks (`samtools quickcheck`), index generation, `flagstat`, `stats`, and target coverage calculation (`samtools bedcov`). |
| `QUALIMAP` | `qualimap` | Generates comprehensive exome quality metrics, GC bias distributions, and insert size profiles. |
| `MULTIQC` | `multiqc` | Aggregates all QC metrics into an interactive, standalone HTML quality report. |
| `INFER_SEX_AND_PLOIDY` | `python` | Infers chromosomal biological sex and allosomal ploidy from chrX and chrY target coverage ratios. |
| `ESTIMATE_CONTAMINATION` | `verifybamid` | Quantifies sample cross-contamination using target allele frequencies. |

---

## 3. Recommended Thresholds

| Metric | Modern High-Density Exomes | Sanger Capture Exomes | Action on Failure |
| :--- | :---: | :---: | :--- |
| **Breadth at >= 10x** | $\ge 90.0\%$ | $\ge 75.0\%$ | Exclude sample or flag as low-confidence |
| **Cohort Correlation ($r$)** | $\ge 0.90$ | $\ge 0.85$ | Exclude from dynamic reference sets |
| **Mean Target Depth** | $\ge 30.0\times$ | $\ge 30.0\times$ | Warning if low depth |
| **Duplicate Read Rate** | $\le 20.0\%$ | $\le 20.0\%$ | Check library complexity |

---

## 4. Output Files
* `qc_results/multiqc_report.html`: Interactive web-based quality summary.
* `qc_results/samtools_metrics/<SAMPLE_ID>.coverage_summary`: Tabular target depth and breadth summary.
* `qc_results/sex_and_ploidy_report.tsv`: Inferred biological sex and karyotypic ploidy table.
