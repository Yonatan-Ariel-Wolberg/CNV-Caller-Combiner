# Evaluation Modules (`modules/evaluation/`)

This directory contains two decoupled evaluation modules, each adhering strictly to the **Single Responsibility Principle (SRP)**:

1. **Performance Metric Evaluation** (`performance-metric-evaluation.nf`): Truthset benchmarking against a gold-standard BED.
2. **Quality Metric Evaluation** (`quality-metric-evaluation.nf`): Callset quality analytics, 3-way comparative matrix, and multi-stage attrition funnel.

---

## Methodology

### 1. Performance Metric Evaluation (`performance-metric-evaluation.nf`)
- **Single Responsibility**: Compares candidate CNV callsets against a reference truthset BED to evaluate detection accuracy.
- **Key Metrics Calculated**:
  - **Sensitivity (Recall)**: Proportion of true CNVs detected.
  - **Precision**: Proportion of called CNVs that match the truthset.
  - **$F_\beta$ Score ($\beta=2$)**: Recall-weighted harmonic mean prioritizing sensitivity over precision.
  - **MCC (Matthews Correlation Coefficient)**: Balanced metric accounting for true and false positives and negatives.
  - **Specificity**: True negative rate where negative spaces are defined.
- **Evaluated Pipeline Stages & Strata**:
  - `Raw CNV Callset` (unfiltered caller outputs)
  - `Quality Filtered CNV Callset` (per-caller score and size filters)
  - `Normalised CNV Callset` (score-calibrated into $[0, 1]$ range)
  - `Intra-Sample Merged CNV Callset` (multi-caller consensus per sample)
  - `Breakpoint-Resolved CNV Callset` (refined coordinates using high-coverage callers)
  - `Cohort-Merged CNV Callset` (cross-sample unified callset)
  - `Allele Frequency-Filtered CNV Callset` (filtered by internal cohort frequency)
  - `CatBoost in silico Validated CNV Callset` (ML-scored and probability-filtered)
  - `ACMG Classified Callsets` (`ACMG_4_5` pathogenic/likely pathogenic, `ACMG_5`, `ACMG_4`, `ACMG_3`, `ACMG_2`, `ACMG_1`, `ACMG_NA`)

### 2. Quality Metric Evaluation (`quality-metric-evaluation.nf`)
- **Single Responsibility**: Extracts intrinsic callset quality metrics and monitors variant filtering attrition across stages without requiring a gold-standard truthset.
- **Variant Attrition Funnel**: Tracks variant retention, drop counts, and percentage remaining from raw inputs down to validated pathogenic variants.
- **3-Way Comparative Matrix**: Evaluates and contrasts every quality metric across:
  1. **Total CNVs**
  2. **Deletions (DEL)**
  3. **Duplications (DUP)**
- **Quality Metrics Extracted per Stage and Subset**:
  - **Call Count & Percentage**: Total variants and relative proportion of stage.
  - **Genomic Footprint**: Total non-redundant megabases (Mb) spanned.
  - **Length Distributions**: Median length (bp), mean length (bp), IQR (bp), min / max length (bp).
  - **Probe Support**: Median probes/call, mean probes/call, single-probe call %, multi-probe ($\ge 3$) call %.
  - **Confidence Scores**: Median QUAL score, mean QUAL score.
  - **Sample Burden**: Average and median calls per sample across the cohort.
- **Output Reports**:
  - `all_stages_quality_report.md`: Markdown summary with 3 clear tables (Total, DEL, DUP).
  - `all_stages_quality_matrix.tsv`: 17-column machine-readable comparative matrix.

---

## Processes in this module

### Performance Evaluation (`performance-metric-evaluation.nf`)
- `EVALUATE_CALLER_PERFORMANCE`: Evaluates sample-level or cohort-level VCF against truthset BED.
- `EVALUATE_ACMG_CSV`: Evaluates ACMG-stratified CSV callsets against truthset BED.
- `EVALUATE_SCORE_SEPARATED`: Evaluates performance stratified by quality score bins.
- `PLOT_EVALUATION_METRICS`: Generates multi-caller ROC and precision-recall curves.

### Quality Analytics (`quality-metric-evaluation.nf`)
- `ANALYZE_CNV_QUALITY_METRICS`: Extracts per-sample quality metrics.
- `ASSESS_CALLSET_QUALITY_METRICS`: Computes the full 3-way Total vs DEL vs DUP metrics for a stage.
- `COMPILE_ALL_STAGES_QUALITY_REPORT`: Compiles multi-stage summary tables and markdown reports.
- `GENERATE_ATTRITION_FUNNEL`: Generates the step-by-step attrition funnel table and chart.

---

## Required parameters
`workflow`: `evaluate`
- `--truth_cnv_bed`: Path to the benchmark truthset BED file (for performance evaluation).
- `--outdir`: Target output directory for evaluation reports and charts.

