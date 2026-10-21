# Nextflow DSL2 Modules (`modules/`)

This directory houses the Nextflow DSL2 process and sub-workflow definitions. Each module encapsulates execution details (containers, environment scripts, and command-line inputs/outputs) for a specific phase of the pipeline.

---

## Complete Inventory of Modules in `modules/`

The pipeline is organized into modular Nextflow DSL2 process definitions under `modules/`:

### 1. CNV & Structural Variant Callers (`callers/`)

- **[canoes.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/canoes.nf)**: Negative Binomial HMM exome CNV caller with GC and length covariates.
- **[clamms.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/clamms.nf)**: Mappability-aware read-depth mixture model caller.
- **[cnvkit.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/cnvkit.nf)**: On-target & antitarget circular binary segmentation (CBS) caller.
- **[exomedepth.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/exomedepth.nf)**: Beta-Binomial caller with dynamic reference cohort selection.
- **[gatk-gcnv.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/gatk-gcnv.nf)**: Hierarchical Bayesian MCMC caller (contig ploidy + GermlineCNVCaller).
- **[icav2-dragen.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/icav2-dragen.nf)**: Illumina Connected Analytics v2 DRAGEN Germline Enrichment caller.
- **[svaba.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/svaba.nf)**: Local de novo string-graph assembly & split-read SV caller.
- **[xhmm.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/xhmm.nf)**: PCA SVD matrix-denoising and Viterbi HMM CNV caller.

### 2. Machine Learning Prioritization (`ml/`)

- **[feature-extraction.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/feature-extraction.nf)**: Extracts 20+ biophysical, log2 ratio, and concordance features per candidate.
- **[cohort-database.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/cohort-database.nf)**: Compiles structured SQLite/DuckDB cohort feature databases.
- **[catboost-train.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/catboost-train.nf)**: Trains CatBoost / XGBoost classifiers and generates SHAP interpretability profiles.
- **[catboost-predict.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/catboost-predict.nf)**: Evaluates and filters candidate CNVs using trained models and Bayesian-tuned thresholds.

### 3. Filtering & Stratification (`filters/`)

- **[quality-filters.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/filters/quality-filters.nf)**: Caller-specific Phred, Bayes Factor, Q_SOME, and probe-count gating (`balanced`, `high_precision`, `max_sensitivity`).
- **[allele-frequency-filters.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/filters/allele-frequency-filters.nf)**: Internal cohort allele frequency filtering ($\le 5\%$) eliminating common technical artifacts.

### 4. SV Merging & Multi-Caller Consensus (`sv-mergers/`)

- **[survivor.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/survivor.nf)**: Reciprocal overlap clustering and intra-/cross-caller union merging.
- **[truvari.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/truvari.nf)**: Sequence- and breakpoint-aware variant collapsing and consensus evaluation.
- **[common-merger.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/common-merger.nf)**: Common merger abstraction for dual Truvari / SURVIVOR routing.

### 5. Sample Quality Assessment (`sample-qc/`)

- **[target-aware-exome-sequencing-quality-assess.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sample-qc/target-aware-exome-sequencing-quality-assess.nf)**: Target-aware depth, mean coverage, Qualimap metrics, sex/contig ploidy inference, and MultiQC aggregation.

### 6. Clinical & Functional Annotation (`annotation/`)

- **[annotation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/annotation/annotation.nf)**: Automated AnnotSV annotation and ACMG/ClinGen 5-class pathogenicity tiering.
- **[functional-annotation-scoring.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/annotation/functional-annotation-scoring.nf)**: StrVCTVRE pathogenicity scoring and ranking.
- **[african-baseline-annotation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/annotation/african-baseline-annotation.nf)**: Baseline population allele frequency contrast.
- **[phen2gene.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/annotation/phen2gene.nf)**: Phenotype-driven HPO gene prioritization.

### 7. Diagnostic Visualization (`visualization/`)

- **[seenv.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/visualization/seenv.nf)**: Cohort-contrast infographic cards and exonic read-depth profiles via SeeNV.
- **[igv.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/visualization/igv.nf)**: Headless, standalone zero-data-leakage IGV HTML reports.

### 8. Metric Evaluation & Benchmarking (`evaluation/`)

- **[performance-metric-evaluation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/evaluation/performance-metric-evaluation.nf)**: Sensitivity, Precision, $F_2$, MCC, and specificity against Sanger truth sets across size/type strata.
- **[quality-metric-evaluation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/evaluation/quality-metric-evaluation.nf)**: DEL/DUP ratios, length distributions, and multi-stage variant attrition funnels.

### 9. Format Conversion & Breakpoint Refinement

- **[format-conversions.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/conversions/format-conversions.nf)**: 40-point forensic CRAM-to-BAM/BAM-to-CRAM conversions and target slicing.
- **[breakpoint-refinement.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/refinement/breakpoint-refinement.nf)**: Exon-bounded read-depth consensus reconciliation with split-read assembly breakends.
- **[normalise-quality-scores.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/normalise/normalise-quality-scores.nf)**: Standardizes quality scores and formats across callers.
- **[bgzip-sort-index.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/common/bgzip-sort-index.nf)**: Coordinate sorting, INFO/TOOL annotation, bgzip compression, and tabix indexing.
