# CNVkit (Copy Number Variation Toolkit)

## 1. Methodology & Algorithm
**CNVkit** is an advanced copy-number caller designed specifically for targeted hybridization capture sequencing. Unlike callers that solely evaluate on-target captured exons, CNVkit incorporates both **on-target** and **off-target** reads:
* **Antitarget Generation**: Constructs off-target "antitarget" bins across non-captured genomic space (typically 100 kb bins), capturing reads from non-specific binding and surrounding regions.
* **Dual Read-Depth Accounting**: Exploits off-target coverage to establish continuous copy-number profiles across large intergenic regions, resolving broad megabase-scale CNVs and whole-chromosome aneuploidies.
* **Bias Correction**: Simultaneously corrects for GC content, sequence repeat content, target bin size, and genomic mappability using median-polish normalization against a pooled cohort reference.
* **Segmentation**: Applies Circular Binary Segmentation (CBS) or Hidden Markov Models to partition normalized log2 depth ratios into distinct discrete copy-number segments.
* **Call Genotyping**: Converts segment log2 ratios to integer copy numbers (CN = 0, 1, 2, 3, 4+) based on sample purity and ploidy expectations.

---

## 2. Nextflow Implementation & Processes
The CNVkit module is encapsulated in `modules/callers/modules-cnvkit.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `CNVKIT_ACCESS` | `cnvkit` | Computes accessible genomic regions from reference FASTA, excluding sequencing assembly gaps (telomeres/centromeres). |
| `CNVKIT_TARGET_ANTITARGET` | `cnvkit` | Generates on-target and off-target antitarget interval files from target BED and accessible regions. |
| `CNVKIT_COVERAGE` | `cnvkit` | Measures on-target (`.targetcoverage.cnn`) and off-target (`.antitargetcoverage.cnn`) read depths per sample. |
| `CNVKIT_REFERENCE` | `cnvkit` | Builds a pooled cohort reference baseline from control samples or flat-reference baseline. |
| `CNVKIT_FIX` | `cnvkit` | Normalizes sample coverage ratios against reference baseline, correcting for GC and mappability biases (`.cnr`). |
| `CNVKIT_SEGMENT` | `cnvkit` | Executes Circular Binary Segmentation (CBS) across log2 ratios to determine copy-number segments (`.cns`). |
| `CNVKIT_CALL` | `cnvkit` | Assigns absolute integer copy numbers and generates VCF formatted output (`.vcf`). |
| `BGZIP_SORT_INDEX_VCF` | `bcftools` | Compresses, sorts, and indexes per-sample VCF files. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA file and `.fai` index.

---

## 4. Key Parameters & Quality Thresholds
* `params.cnvkit_method`: Segmentation method (default: `'cbs'`, option: `'hmm'`).
* `params.filter_cnvkit_qscore`: Quality score threshold (default: `>= 30`).
* `params.filter_cnvkit_probes`: Minimum number of supporting intervals (default: `>= 3`).
* `params.filter_cnvkit_min_size`: Minimum length in base pairs (default: `1000 bp`).

---

## 5. Output Files
* `out_CNVKIT/raw/<SAMPLE_ID>.cnvkit.raw.vcf.gz`: Unfiltered VCF containing all candidate segments.
* `out_CNVKIT/filtered/<SAMPLE_ID>.cnvkit.filtered.vcf.gz`: Quality-filtered and normalized callset.
* `out_CNVKIT/segments/<SAMPLE_ID>.cns`: Detailed segment log2 ratios and integer copy numbers.
