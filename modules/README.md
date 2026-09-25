# Nextflow DSL2 Modules (`modules/`)

This directory houses the Nextflow DSL2 process and sub-workflow definitions. Each module encapsulates execution details (containers, environment scripts, and command-line inputs/outputs) for a specific phase of the pipeline.

---

## 📁 Directory Structure

### 1. Caller Pipelines (`modules/callers/`)
Contains process modules for individual CNV calling tools:
- **[modules-canoes.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-canoes.nf)**: Execution of CANOES.
- **[modules-clamms.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-clamms.nf)**: Execution of CLAMMS.
- **[modules-cnvkit.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-cnvkit.nf)**: Execution of CNVkit.
- **[modules-gatk-gcnv.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-gatk-gcnv.nf)**: Execution of GATK gCNV.
- **[modules-icav2-dragen.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-icav2-dragen.nf)**: Execution of ICAv2 cloud-based DRAGEN CNV calling.
- **[modules-exomedepth.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-exomedepth.nf)**: Execution of ExomeDepth.
- **[modules-xhmm.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-xhmm.nf)**: Execution of XHMM.

### 2. Common Operations (`modules/common/`)
General utility workflows:
- **[bgzip-sort-index.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/common/bgzip-sort-index.nf)**: Standard per-sample VCF coordinate-sorting, tool annotation (`INFO/TOOL`), bgzip compression, and tabix indexing.
- **[modules-common.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/common/modules-common.nf)**: Backwards-compatible common utility re-exports.
- **[modules-mosdepth-target.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/common/modules-mosdepth-target.nf)**: High-performance target-aware depth estimation using Mosdepth.

### 3. Metric Evaluation (`modules/evaluation/`)
- **[performance-metric-evaluation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/evaluation/performance-metric-evaluation.nf)**: Evaluates CNV callsets at different stages against gold-standard truth BED in terms of Sensitivity, Precision, F-beta (beta=2), MCC, and Specificity, comparing performance across all stages and strata.
- **[quality-metric-evaluation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/evaluation/quality-metric-evaluation.nf)**: Analyzes callset quality metrics (DEL/DUP ratios, length distributions, footprint Mb, median probes/call, sample calls), generates the variant attrition funnel, and compiles multi-stage summary tables.

### 4. Machine Learning & In Silico Validation (`modules/ml/`)
- **[modules-feature-extraction.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/modules-feature-extraction.nf)**: Runs feature engineering steps on merged candidate sites.
- **[modules-train.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/modules-train.nf)**: Trains classifiers and exports SHAP interpretability profiles.
- **[modules-predict.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/modules-predict.nf)**: Applies trained CatBoost models to score and filter candidate CNVs.
- **[modules-cohort-database.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/ml/modules-cohort-database.nf)**: Compiles structured SQLite/DuckDB cohort feature databases.

### 5. Filters (`modules/filters/`)
- **[quality-filters.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/filters/quality-filters.nf)**: Quality and length filters for individual caller callsets, dynamically controlled by `--filter_profile` (`balanced`, `high_precision`, `max_sensitivity`) or individual CLI overrides.
- **[allele-frequency-filters.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/filters/allele-frequency-filters.nf)**: Filters cohort callsets by internal dataset allele frequency (AF <= 5%), strictly focused internally on cohort samples without external databases (ClinVar, gnomAD, etc.).

### 6. Annotation & ACMG Classification (`modules/annotation/`)
- **[annotation.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/annotation/annotation.nf)**: AnnotSV database setup, external annotation merging, and ACMG pathogenicity classification/stratification.

### 7. Normalization (`modules/normalise/`)
- **[normalise-quality.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/normalise/normalise-quality.nf)**: Normalises CNV caller quality scores to a common standardized scale downstream of quality filtering.
- **[modules-normalise.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/normalise/modules-normalise.nf)**: Standalone wrapper for normalising pre-existing caller VCF channels.

### 8. Sample Quality Control (`modules/sample-qc/`)
- **[target-aware-exome-sequencing-quality-assess.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sample-qc/target-aware-exome-sequencing-quality-assess.nf)**: Executes exome coverage assessment to flags samples with insufficient depth or skewed metrics.

### 9. SV Merging & Consensus (`modules/sv-mergers/`)
Orchestrates sample-level and cohort-level callset combination:
- **[modules-survivor.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/modules-survivor.nf)**: Runs SURVIVOR for genomic overlap merging.
- **[modules-truvari.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/modules-truvari.nf)**: Runs Truvari for VCF-based matching and consensus.
- **[modules-common-merger.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/modules-common-merger.nf)** / **[modules-truvari-prefilter.nf](file:///home/ywolberg/CNV-Caller-Combiner/modules/sv-mergers/modules-truvari-prefilter.nf)**: Helpers for pre-filtering VCF variants prior to combination.
