# GATK gCNV (Germline CNV Caller)

## 1. Methodology & Algorithm
The **Broad Institute GATK GermlineCNVCaller (gCNV)** is a state-of-the-art probabilistic framework for discovering copy-number variation from targeted sequencing depth:
* **Hierarchical Bayesian Model**: Models read-depth observations as negative-binomial distributed variables conditioned on latent interval copy-number states, sample-specific sequencing depth, and continuous interval-specific bias factors (GC content, mappability, probe effects).
* **Markov Chain Monte Carlo (MCMC) & Variational Inference**: Infers model parameters simultaneously across samples and intervals, learning systematic noise components without requiring external controls.
* **Cohort Contig Ploidy Modeling**: Explicitly models allosome (chrX/chrY) and autosome copy-number ploidy per sample before interval calling using `DetermineGermlineContigPloidy`.
* **Posterior Genotype Probabilities**: Assigns precise posterior probabilities to each copy number state, generating quality metrics (`QS`, `QSS`, and Phred-scaled likelihoods).

---

## 2. Nextflow Implementation & Processes
The GATK gCNV pipeline module is encapsulated in `modules/callers/modules-gatk-gcnv.nf` and contains the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `GATK_PREPROCESS_INTERVALS` | `gatk` | Preprocesses and pads target BED intervals into Picard-style `.interval_list` format. |
| `GATK_COLLECT_READ_COUNTS` | `gatk` | Collects per-interval read counts from BAM/CRAM alignments in `HDF5` format. |
| `GATK_DETERMINE_PLOIDY` | `gatk` | Determines sample-level karyotypic sex and baseline contig ploidy distributions across the cohort. |
| `GATK_SCATTER_INTERVALS` | `gatk` | Scatters interval lists across multiple shards for parallel Bayesian model fitting. |
| `GATK_GERMLINE_CNV_CALLER` | `gatk` | Fits the hierarchical Bayesian model and infers copy-number states per shard. |
| `GATK_POSTPROCESS_CALLS` | `gatk` | Merges sharded calls into unified per-sample VCFs with genotype likelihoods. |
| `FORMAT_GCNV_VCF` | `python` | Standardizes GATK gCNV output attributes for consensus merging. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM files (`.bam` + `.bai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA file with `.fai` and `.dict` files.
* **Ploidy Priors Table**: Baseline contig ploidy prior specifications for autosomes and allosomes.

---

## 4. Key Parameters & Quality Thresholds
* `params.filter_gcnv_qscore`: Quality score threshold (default: `>= 30`).
* `params.filter_gcnv_probes`: Minimum number of supporting intervals (default: `>= 3`).
* `params.filter_gcnv_min_size`: Minimum event length in base pairs (default: `1000 bp`).

---

## 5. Output Files
* `out_GCNV/raw/<SAMPLE_ID>.gcnv.raw.vcf.gz`: Unfiltered VCF containing all candidate segments.
* `out_GCNV/filtered/<SAMPLE_ID>.gcnv.filtered.vcf.gz`: Quality-filtered and normalized callset.
