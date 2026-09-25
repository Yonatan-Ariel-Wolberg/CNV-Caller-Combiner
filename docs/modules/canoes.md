# CANOES (Calculated Armature for Novel Variant Occurrence in Exome Sequencing)

## 1. Methodology & Algorithm
**CANOES** is a read-depth based germline Copy Number Variation (CNV) caller specifically engineered for whole-exome sequencing (WES) data. Because WES targets discrete captured regions (exons) rather than uniform whole genomes, read counts exhibit severe target-to-target variance driven by GC content fluctuations, probe hybridization affinities, and regional mappability.

CANOES addresses this by modeling the read counts of each targeted interval across a cohort of samples using a **Negative Binomial distribution** coupled with a **Hidden Markov Model (HMM)**:
* **Negative Binomial Modeling**: Accounts for the substantial over-dispersion inherent to next-generation sequencing count data (where variance exceeds the mean), which standard Poisson models underestimate.
* **GC-Content & Depth Covariates**: Explicitly normalizes per-interval read counts against reference GC percentages and whole-sample global sequencing depth.
* **Viterbi Decoding**: Leverages dynamic programming across genomic coordinates to identify contiguous runs of exons exhibiting copy loss (deletion: CN = 0 or 1) or copy gain (duplication: CN >= 3).
* **Likelihood Calculations**: Assigns quality scores (MLCNV score) proportional to the log-likelihood difference between the CNV state and the normal diploid state.

---

## 2. Nextflow Implementation & Processes
The CANOES pipeline module is encapsulated in `modules/callers/modules-canoes.nf` and includes the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `CALC_GC_CANOES` | `gatk` | Computes interval GC percentages per autosome (`chr1`-`chr22`) using `gatk AnnotateIntervals` on the reference FASTA. |
| `BATCH_BAMS` | `canoes_R` | Splits cohort BAM lists into fixed-size batches (default: 100 BAMs) to optimize memory usage and prevent file descriptor limits. |
| `CANOES_READS` | `bedtools` | Quantifies interval read depths in parallel batches across target BED intervals using `bedtools multicov`. |
| `COMBINE_READS` | `samtools` | Concatenates batch-wise read depth tables into unified per-chromosome count matrices. |
| `RUN_CANOES` | `canoes_R` | Executes the core CANOES R statistical package across samples and autosomes, fitting the negative binomial models and outputting raw CNV segments. |
| `FORMAT_CANOES_VCF` | `python` | Converts tabular CANOES calls into standardized VCF 4.2 format, computing `SVLEN`, `END`, and Phred-scaled confidence scores. |
| `BGZIP_SORT_INDEX_VCF` | `bcftools` | Sorts, bgzips, and indexes the resulting per-sample VCFs with `tabix`. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons (3-column or 4-column BED format).
* **Reference Genome**: GRCh38/hg38 FASTA file and companion `.fai` index.
* **Cohort Size**: Recommended minimum cohort size of $\ge 15$ samples to provide sufficient statistical power for the negative binomial dispersion parameter estimation.

---

## 4. Key Parameters & Quality Thresholds
* `params.canoes_batch_size`: Batch size for parallel `bedtools multicov` execution (default: `100`).
* `params.filter_canoes_qscore`: Phred-scaled quality score cutoff (default: `>= 30` in balanced mode).
* `params.filter_canoes_probes`: Minimum number of targeted exons required to support a CNV event (default: `>= 3` probes).

---

## 5. Output Files
* `out_CANOES/raw/<SAMPLE_ID>.canoes.raw.vcf.gz`: Unfiltered VCF containing all candidate deletions and duplications.
* `out_CANOES/filtered/<SAMPLE_ID>.canoes.filtered.vcf.gz`: Quality-gated VCF filtered by caller score and probe count.
