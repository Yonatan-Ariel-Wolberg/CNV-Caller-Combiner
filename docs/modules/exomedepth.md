# ExomeDepth (Beta-Binomial Exome CNV Calling)

## 1. Methodology & Algorithm
**ExomeDepth** is a widely acclaimed statistical method specifically engineered to call germline copy number variants from targeted exome sequencing data:
* **Beta-Binomial Read-Count Distribution**: Explicitly accounts for technical variation in read counts across targeted intervals. While Poisson models assume variance equals the mean, ExomeDepth models the over-dispersion inherent to targeted sequencing by using a beta-binomial distribution where the probability of a read mapping to an exon varies across samples.
* **Dynamic Reference Cohort Selection**: Rather than using a static reference panel or averaging all samples, ExomeDepth uses a forward-selection algorithm to build an optimized reference set for each individual sample from the remaining cohort samples. It selects samples that have the highest pairwise correlation (typically requiring $r \ge 0.97$), ensuring the reference closely mirrors the test sample's specific capture efficiency and technical batch profile.
* **Hidden Markov Model (HMM)**: Uses an HMM across exons along each chromosome to distinguish true multi-exon copy-number events from single-exon technical fluctuations.
* **Bayes Factor (BF) Scoring**: Calculates a Bayes Factor comparing the likelihood of a CNV model against the normal diploid model. A high Bayes Factor (typically $\ge 10.0$ or higher) denotes strong statistical support for a true variant.

---

## 2. Nextflow Implementation & Processes
The ExomeDepth module is encapsulated in `modules/callers/modules-exomedepth.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `PREPARE_EXOMEDEPTH_TARGETS` | `bedtools` | Generates a 4-column BED file (chrom, start, end, name) formatted for ExomeDepth interval accounting. |
| `BATCH_BAMS_EXOMEDEPTH` | `canoes_R` | Splits cohort BAM lists into manageable batches (default: 50 BAMs) to distribute computation and prevent file handle exhaustion. |
| `EXOMEDEPTH_READS_BATCH` | `bedtools` | Quantifies target exon read counts across each batch in parallel using `bedtools multicov`. |
| `COMBINE_EXOMEDEPTH_READS` | `python` | Merges per-batch read count tables into a single unified cohort-wide read count matrix. |
| `RUN_EXOMEDEPTH_CALLING` | `canoes_R` | Executes the ExomeDepth R package: performs dynamic reference selection, fits beta-binomial distributions, and decodes CNV calls via HMM. |
| `FORMAT_EXOMEDEPTH_VCF` | `python` | Converts tabular ExomeDepth calls into standard VCF 4.2 format, populating `SVLEN`, `END`, `BF` (Bayes Factor), and Phred-scaled `QUAL`. |
| `BGZIP_SORT_INDEX_VCF` | `bcftools` | Sorts, bgzips, and creates `.tbi` indices for downstream merging. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA and `.fai` index.
* **Cohort Size**: Requires at least 10–15 samples to enable the dynamic reference selection algorithm to find suitably correlated reference matches ($r \ge 0.97$).

---

## 4. Key Parameters & Quality Thresholds
* `params.exomedepth_batch_size`: Number of BAMs processed per `bedtools multicov` task (default: `50`).
* `params.exomedepth_min_bf`: Minimum Bayes Factor required for a call (default: `10.0`).
* `params.exomedepth_min_correlation`: Minimum correlation cutoff when selecting reference samples (default: `0.97`).
* `params.filter_exomedepth_qscore`: Phred-scaled quality score threshold (default: `>= 30`).
* `params.filter_exomedepth_probes`: Minimum number of supporting exons (default: `>= 3`).

---

## 5. Output Files
* `out_EXOMEDEPTH/raw/<SAMPLE_ID>.exomedepth.raw.vcf.gz`: Raw VCF calls with Bayes Factor annotations.
* `out_EXOMEDEPTH/filtered/<SAMPLE_ID>.exomedepth.filtered.vcf.gz`: Filtered VCF calls passing quality thresholds.
