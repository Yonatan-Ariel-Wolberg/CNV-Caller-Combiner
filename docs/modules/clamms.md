# CLAMMS (Coverage Loss Automation for Microarray and Mass Sequencing)

## 1. Methodology & Algorithm
**CLAMMS** is an exome-optimized germline CNV calling algorithm designed to deliver high sensitivity and low false-discovery rates by leveraging pre-computed regional genomic features:
* **Mappability Modeling**: Integrates 75 bp, 100 bp, or 150 bp read-length k-mer mappability tracks to account for repetitive genomic architecture, pseudogenes, and segmental duplications.
* **Specialized Reference Models**: Constructs specialized baseline models for custom assay designs, adjusting target interval depths using Picard insert size metrics and regional GC composition.
* **Window-Level Depth Normalization**: Normalizes coverage per interval against global coverage distributions, fitting mixture models to detect fractional copy number states.
* **Multi-Allelic Resolution**: Resolves diploid and non-diploid copy numbers (including zero-copy deletions, single-copy deletions, normal diploid states, duplications, and high-copy amplifications).

---

## 2. Nextflow Implementation & Processes
The CLAMMS pipeline module is encapsulated in `modules/callers/modules-clamms.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `PREPARE_CLAMMS_BED` | `bedtools` | Generates standardized interval BED coordinates with annotated GC percentage and mappability scores. |
| `CLAMMS_COVERAGE` | `clamms` | Computes target coverage per sample using BAM/CRAM alignments and specialized interval definitions. |
| `CLAMMS_FIT_MODELS` | `clamms` | Fits statistical baseline coverage models across cohort samples, calculating expected per-window distributions. |
| `CLAMMS_CALL_CNV` | `clamms` | Evaluates coverage deviations against fitted reference models, calling segment-level deletions and duplications. |
| `FORMAT_CLAMMS_VCF` | `python` | Translates raw tabular CLAMMS CNV calls into standard VCF 4.2 format with `SVLEN`, `END`, and log-odds confidence metrics. |
| `NORM_CLAMMS_FULL` | `python` | Standardizes raw CLAMMS quality scores to Phred-scale (0–99). |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA file with `.fai` index.
* **Mappability Tracks**: Pre-computed BigWig/BED mappability files matching sequencing read length (e.g. 100 bp mappability).
* **Insert Size Metrics**: Picard `CollectInsertSizeMetrics` tables per sample.

---

## 4. Key Parameters & Quality Thresholds
* `params.filter_clamms_qscore`: Quality score threshold (default: `>= 30` in balanced mode).
* `params.filter_clamms_probes`: Minimum number of supporting intervals (default: `>= 3`).
* `params.filter_clamms_min_size`: Minimum CNV event length in base pairs (default: `1000 bp`).

---

## 5. Output Files
* `out_CLAMMS/raw/<SAMPLE_ID>.clamms.raw.vcf.gz`: Raw output VCF containing candidate calls.
* `out_CLAMMS/filtered/<SAMPLE_ID>.clamms.filtered.vcf.gz`: Quality-filtered and normalized callset.
