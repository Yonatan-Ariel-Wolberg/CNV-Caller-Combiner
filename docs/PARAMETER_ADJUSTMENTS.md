# Benchmark Parameter Adjustments & Methodological Justifications

This document details all parameter modifications, threshold tunings, heuristic adjustments, and engineering adaptations applied to third-party tools within the **CNV-Caller-Combiner** pipeline. Each section contrasts the tool's standalone or default behavior with the customized pipeline settings and articulates the biological and computational rationale for each change.

---

## 1. Depth-Based Germline CNV Callers

### 1.1 CANOES
*Reference: Backenroth et al., Nucleic Acids Res 2014*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Mapping Quality (`-q`)** | `bedtools multicov -q 0` | `-q 20` | Filters ambiguously mapped reads and low-quality multimappers before regression modeling. |
| **BAM Batching** | Monolithic (all BAMs at once) | `canoes_batch_size = 100` | Prevents operating system `Too many open files` file-descriptor exhaustion during large cohort runs. |
| **Noise Outlier Gating** | Retains all sample calls | Samples with $\ge 50$ calls/chr segregated to `xcnv.fail` | Quarantines hyper-noisy or degraded samples, preventing them from corrupting cohort-wide frequency calculations. |
| **Quality Score (`Q_SOME`)** | Recommended $\ge 80$ | • **Balanced**: `Q_SOME >= 65`<br>• **High Precision**: `Q_SOME >= 80`<br>• **Max Sensitivity**: `Q_SOME >= 50` | Calibrates sensitivity/specificity tradeoff across benchmark arms. |
| **Minimum Event Size** | Unconstrained | • **Balanced**: $500\text{ bp}$<br>• **High Precision**: $1,000\text{ bp}$<br>• **Max Sensitivity**: $0\text{ bp}$ | Prunes sub-microscopic single-probe noise in precision-oriented arms. |
| **Normalized Score Formulation** | None (native Phred) | $\text{QUAL}_{\text{norm}} = \min(1000.0, 100.0 \times (\frac{Q\_SOME}{80.0}))$ | Standardizes caller scores to a uniform [0, 1000] scale for multi-caller consensus. |

---

### 1.2 CLAMMS
*Reference: Packer et al., Bioinformatics 2016*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Coverage Quality (`-Q`)** | `samtools bedcov -Q 0` | `-Q 30` | Ensures high-confidence mapping quality over target exon windows. |
| **Fragment Insert Size** | Variable / uncalibrated | `INSERT_SIZE=200` | Accurately models paired-end exome library fragment length during window generation. |
| **Mappability Window** | Unconstrained | 75-mer mappability track (`mappability_sorted.bed`) | Excludes non-unique and repetitive genomic windows prone to false-positive depth fluctuations. |
| **Likelihood Thresholds (`Q_SOME`, `Q_EXACT`)** | `Q_SOME >= 50`<br>`DEL: Q_EXACT >= 0.0`<br>`DUP: Q_EXACT >= -0.5` | • **Balanced**: `Q_SOME >= 50, DEL >= 0.0, DUP >= -0.5`<br>• **High Precision**: `Q_SOME >= 100, DEL >= 0.5, DUP >= -0.2`<br>• **Max Sensitivity**: `Q_SOME >= 30, DEL >= -1.0, DUP >= -1.0` | Enables distinct, strictly monotonic benchmark evaluation across profiles while accommodating duplication signal-to-noise dynamics. |
| **Sex-Aware Baseline** | Autosomal default | Dynamic `--sex` lookup from cohort metadata | Configures expected baseline copy number on allosomal contigs (chrX/chrY). |
| **Normalized Score Formulation** | None (native log-likelihood) | $\text{QUAL}_{\text{norm}} = \min(1000.0, 100.0 \times (\frac{Q\_SOME}{50.0}))$ | Projects CLAMMS confidence onto universal [0, 1000] scale. |

---

### 1.3 XHMM
*Reference: Fromer et al., Am J Hum Genet 2012*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Coverage Batching** | Monolithic GATK run | `xhmm_batch_size = 50` | Mitigates JVM heap exhaustion and CPU thrashing during multi-sample GATK DepthOfCoverage. |
| **Mapping Quality (`-mmq`)** | `-mmq 0` | `-mmq 20 --count-type COUNT_READS` | Enforces MAPQ $\ge 20$ and counts discrete read alignments rather than fragmented bases. |
| **Extreme GC Pruning** | Manual inspection | Automated filter excluding targets with $\text{GC} < 0.10$ or $\text{GC} > 0.90$ | Prevents extreme GC-biased exons from destabilizing singular value decomposition (SVD) during PCA. |
| **PCA Variance Factor** | Variable | `PVE_mean_factor 0.70` | Strips out systematic batch and capture effects while preserving authentic biological copy-number variance. |
| **Z-Score Outlier Bound** | Unconstrained | `--maxSdTargetRD 30` | Excludes hyper-variable target intervals where read-depth standard deviation exceeds 30. |
| **Quality Thresholds (`SQ`, `NDQ`, `EQ`)** | Recommended `SQ >= 60` | • **Balanced**: `SQ >= 50, NDQ >= 30, EQ >= 30, Size >= 500 bp`<br>• **High Precision**: `SQ >= 60, NDQ >= 60, EQ >= 60, Size >= 1000 bp`<br>• **Max Sensitivity**: `SQ >= 30, NDQ >= 0, EQ >= 0, Size >= 0 bp` | Explores the precision/recall Pareto frontier across filtering regimes. |
| **Normalized Score Formulation** | Native Phred quality | $\text{QUAL}_{\text{norm}} = \min(1000.0, 100.0 \times (\frac{SQ}{60.0}))$ | Scales XHMM Phred quality score against the standard clinical threshold (60.0). |

---

### 1.4 CNVkit
*Reference: Talevich et al., PLOS Comput Biol 2016*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Autobin Dynamic Sizing** | Static 100 kb bins | Dynamic target/antitarget estimation based on median BAM coverage | Automatically adapts bin widths to exome capture kit density and sequencing depth. |
| **Parallel Coverage** | Single-threaded | Multi-threaded (`-p 2` for target and antitarget concurrently) | Maximizes core utilization and halves per-sample coverage collection latency. |
| **Segment Filtering & Probe Gating** | Retains all segments (including 1-probe calls) | • **Balanced**: $\ge 3\text{ probes}, \ge 5\text{ kb}$<br>• **High Precision**: $\ge 5\text{ probes}, \ge 10\text{ kb}$<br>• **Max Sensitivity**: $\ge 2\text{ probes}, \ge 1\text{ kb}$ | Excludes isolated single-exon false positives that plague unconstrained exome segmentation. |
| **Quality & Size Bounds** | Unbounded | `FMT/GQ >= 1 || FMT/CNQ >= 1`, max length capped at $5\text{ Mb}$ | Purges unsegmented telomeric and centromeric noise artifacts. |
| **Normalized Score Formulation** | Log2 ratio & CNQ | $\text{QUAL}_{\text{norm}} = \min(1000.0, 100.0 \times (\frac{CNQ}{15.0}) + \min(100.0, |\text{log2}| \times 50.0))$ | Rewards high copy-number likelihood (CNQ) and large absolute log2-ratio effect sizes. |

---

### 1.5 GATK gCNV (GermlineCNVCaller)
*Reference: Broad Institute Germline CNV Best Practices*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Exome Target Padding** | 0 bp padding | 250 bp padding (`--padding 250`, `--bin-length 0`) | Captures flanking read depth around exon-intron boundaries while maintaining native probe boundaries. |
| **Coherence Lengths** | `10000.0 bp` (WGS default) | `1000.0 bp` (`--class-coherence-length 1000.0 --cnv-coherence-length 1000.0`) | Prevents the hidden Markov model from oversmoothing focal, single-exon deletions and duplications. |
| **Hyperparameter Tunings** | Standard defaults | `--enable-bias-factors false --interval-psi-scale 1.0E-6 --log-mean-bias-standard-deviation 0.01 --sample-psi-scale 1.0E-6` | Stabilizes variance parameters for targeted capture exomes without over-fitting. |
| **Contig Ploidy Priors** | Uniform priors across contigs | Tailored priors for autosomes, chrX, and chrY in mixed cohorts | Accurately models hemizygous sex chromosomes in males vs. females. |
| **Interval Sharding** | Monolithic | `scatter_count = 5000` (via `IntervalListTools`) | Distributes compute across cluster nodes, reducing runtime from days to hours. |
| **Postprocessing Quality Thresholds** | Unconstrained segments | • **Balanced**: `QS_DEL >= 100, QS_DUP >= 50, single_interval_qs >= 100`<br>• **High Precision**: `QS_DEL >= 400, QS_DUP >= 400, single_interval_qs >= 1000`<br>• **Max Sensitivity**: `QS_DEL >= 50, QS_DUP >= 20, single_interval_qs >= 50` | Enforces stricter gating on single-interval calls ($NP < 2$) to eliminate capture probe noise. |

---

### 1.6 ExomeDepth
*Reference: Plagnol et al., Bioinformatics 2012*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Reference Correlation Cutoff** | Flexible ($r \ge 0.90$) | $r \ge 0.97$ (`min_correlation = 0.97`) | Enforces high technical correlation among reference aggregate exomes to minimize background noise. |
| **Transition Probability** | Variable | Fixed at $10^{-4}$ (`transition_prob = 1e-4`) | Balances sensitivity for focal CNVs against over-segmentation of single events. |
| **Multicov Batching** | Single run | Batched to 50 BAMs per job (`exomedepth_batch_size = 50`) | Prevents file descriptor exhaustion during cohort count matrix aggregation. |
| **Bayes Factor (`BF`) & Size Filtering** | Package default reports all calls | • **Balanced**: `BF >= 10.0, Size >= 500 bp`<br>• **High Precision**: `BF >= 30.0, Size >= 1000 bp`<br>• **Max Sensitivity**: `BF >= 5.0, Size >= 0 bp` | Separates marginal calls from high-confidence clinical candidates. |
| **Normalized Score Formulation** | Bayes Factor (`BF`) | $\text{QUAL}_{\text{norm}} = \min(1000.0, 100.0 \times (\frac{BF}{20.0}))$ | Standardizes Bayes Factor scores against the benchmark threshold ($BF = 20.0$). |

---

### 1.7 SvABA
*Reference: Wala et al., Genome Res 2018*

| Parameter / Configuration | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Target Restriction** | Whole genome assembly | `-k ${target_bed}` | Restricts local de novo string graph assembly strictly to capture targets, cutting compute time from 12+ hrs to ~15 min per sample. |
| **Calling Mode** | Somatic tumor-normal | `--germline` | Configured for constitutional germline structural variant discovery in non-cancer cohorts. |
| **Reference Indexing** | Requires existing BWA index | Automated `.bwt` check & on-demand symlinking | Guarantees required BWA indices are available without manual user intervention. |

---

## 2. Structural Variant Merging & Consensus Engines

### 2.1 Truvari Collapse (Intra-Sample Consensus)
*Reference: English et al., Genome Biology 2022*

| Parameter / Option | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Execution Mode** | Inter-sample comparison | `--intra` | Restricts collapsing to calls within the same individual sample across multiple callers. |
| **Sequence Similarity (`--pctseq`)** | `--pctseq 0.70` (70% identity) | `--pctseq 0` | Read-depth callers output symbolic `<DEL>`/`<DUP>` intervals; sequence comparison is not applicable and would prevent all matches. |
| **Size Similarity (`--pctsize`)** | `--pctsize 0.70` | `--pctsize 0.50` | Allows merging across tools with disparate breakpoint estimation precision (e.g. CBS vs. HMM). |
| **Reciprocal Overlap (`--pctovl`)** | `--pctovl 0.0` | `--pctovl 0.50` | Requires at least 50% coordinate overlap between callers before declaring consensus. |
| **Size Constraints** | Min 50 bp, Max unconstrained | `--sizemin 50 --sizemax 2000000` | Filters sub-50 bp micro-indels and caps upper bound at 2 Mb to exclude gross cytogenetic anomalies. |
| **Clustering Mode** | Standard graph | `--fast-cluster` | Optimizes memory usage and speeds up graph component resolution. |

---

### 2.2 SURVIVOR Merge
*Reference: Jeffares et al., Nat Commun 2017*

| Parameter / Option | Default (Standalone) | Workflow Adjustment | Methodological Rationale |
| :--- | :--- | :--- | :--- |
| **Distance Window (`max_dist`)** | 1,000 bp | • **Intra-Sample**: `1,000 bp`<br>• **Cohort-Level**: `10,000 bp` | 1 kb tightly pairs cross-caller consensus within a sample; 10 kb groups recurrent cohort alleles with staggered probe boundaries. |
| **Strand Requirement (`use_strand`)** | `1` (strand matching required) | `0` (strand matching disabled) | Read-depth CNV callers do not provide strand polarity; enforcing strand matching causes false negatives. |
| **Type Agreement (`use_type`)** | `1` | `1` | Strictly enforces matching SV types (deletions only merge with deletions; duplications with duplications). |
| **Minimum Event Size (`min_sv_size`)** | 0 bp | `30 bp` | Eliminates spurious sub-exon artifacts $< 30$ bp. |
| **Genotype Format Normalization** | Raw caller genotype | Automated haploid-to-diploid conversion (`1` $\rightarrow$ `1/1`) | Prevents SURVIVOR's 3-character fixed-width genotype parser bug from splitting colon fields and corrupting output VCFs. |

---

## 3. Post-Consensus Breakpoint Refinement (`refine_breakpoints.py`)

*Engine: In-house post-consensus structural variant breakpoint reconciler.*

| Parameter / Threshold | Pipeline Setting | Methodological Rationale |
| :--- | :--- | :--- |
| **Maximum Flank Window (`--max-flank`)** | `50,000 bp` | Searches up to 50 kb into unsequenced flanking introns for orthogonal split-read/assembly breakpoints (from SvABA/DRAGEN). |
| **Displacement Delta (`--max-delta`)** | `25,000 bp` | Caps the permissible coordinate shift between a depth-inferred probe boundary and a junction breakpoint at 25 kb. |
| **Focal Threshold (`--focal-threshold`)** | `10,000 bp` | Distinguishes focal events ($\le 10\text{ kb}$) from multi-exon/large CNVs to apply tighter displacement tolerance. |
| **Minimum Focal Delta (`--min-focal-delta`)**| `1,000 bp` | Prevents over-adjusting small focal events with ambiguous alignment junctions. |
| **Reciprocal Overlap (`--min-recip-overlap`)**| `0.20` (20%) | Requires at least 20% mutual overlap to prevent pairing distinct neighboring CNVs within the same locus. |

---

## 4. Cohort QC & Consensus Target BED Generation

### 4.1 Target-Aware Cohort QC
*Reference: Target-Aware QC guidelines for diverse and clinical cohorts.*

| Metric / Threshold | Pipeline Cutoff | Methodological Rationale |
| :--- | :--- | :--- |
| **Breadth of Coverage ($\ge 10\times$)** | $\ge 75.0\%$ | Calibrated to retain clinically valuable samples with high genetic diversity while gating failed capture libraries. |
| **PCR Duplication Rate** | $\le 20.0\%$ | Quarantines low-complexity or over-amplified libraries that produce artificial depth spikes. |
| **Average On-Target Depth** | $\ge 20.0\times$ | Enforces baseline coverage required for statistical power in negative binomial and beta-binomial models. |

### 4.2 Mosdepth Consensus Target BED Generation
*Reference: Pedersen et al., Bioinformatics 2018*

| Parameter / Operation | Pipeline Setting | Methodological Rationale |
| :--- | :--- | :--- |
| **Per-Sample Coverage Threshold** | $\ge 20.0\times$ mean depth per interval | Identifies intervals with sufficient coverage for reliable CNV calling. |
| **Cohort Frequency Threshold** | $\ge 75\%$ passing samples (`pass_percent = 75`) | Excludes kit-specific uncaptured regions while retaining consistent cohort-wide targets. |
| **Gap Bridging (`bedtools merge -d`)** | `1,500 bp` (`-d 1500`) | Bridges unsequenced intronic probe gaps $< 1.5$ kb, preventing artificial split-calling across multi-probe exons. |

---

## 5. Machine Learning & In Silico Validation (CatBoost)

| Parameter / Configuration | Pipeline Setting | Methodological Rationale |
| :--- | :--- | :--- |
| **Classification Model** | CatBoostClassifier | Optimized for tabular datasets containing mixed categorical (caller name, SV type) and numerical features. |
| **Decision Threshold** | `0.50` (`--threshold 0.5`) | Calibrated on clinical truth sets to balance precision against sensitivity for rare pathogenic variants. |
| **Feature Set** | Concordance, Normalized QUAL, Mean QUAL, Best Tool Score, GC Content, Exon Count | Combines caller agreement with underlying genomic architecture features. |

---

## 6. Clinical Annotation & Variant Prioritization (AnnotSV)

*Reference: Geoffroy et al., Bioinformatics 2018; Rigaux et al., Eur J Hum Genet 2021*

| Option / Parameter | Pipeline Setting | Methodological Rationale |
| :--- | :--- | :--- |
| **Annotation Mode** | `-overwrite 0` | Reuses static pre-compiled Human annotation tables, reducing runtime by 80%. |
| **Genome Build** | `GRCh38` (`-genomeBuild GRCh38`) | Aligns with the GRCh38 full analysis set reference genome. |
| **ACMG Classification** | Classes 1 through 5 | Implements ClinGen/ACMG guidelines, stratifying candidates and isolating Class 4 (Likely Pathogenic) and Class 5 (Pathogenic) variants. |
