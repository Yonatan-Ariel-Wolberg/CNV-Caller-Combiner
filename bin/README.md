# Helper Scripts and Executables (`bin/`)

This directory contains the custom scripts, helper tools, plotting utilities, and container recipe definitions invoked by the Nextflow pipeline processes.

---

## 🛠️ Script Categories

### 1. CNV Caller Helper Scripts
- **[CANOES.R](file:///home/ywolberg/CNV-Caller-Combiner/bin/CANOES.R)**: Core CANOES algorithm code.
- **[run_canoes.R](file:///home/ywolberg/CNV-Caller-Combiner/bin/run_canoes.R)**: R script wrapper to execute CANOES calling.
- **[custom_ref_panel.R](file:///home/ywolberg/CNV-Caller-Combiner/bin/custom_ref_panel.R)** / **[custom_ref_panel.sh](file:///home/ywolberg/CNV-Caller-Combiner/bin/custom_ref_panel.sh)**: Compile reference panels for callers.

### 2. Format Converters & Quality Normalizers
- **[canoes_csv_to_vcf.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/canoes_csv_to_vcf.py)**: Converts raw CANOES CSV outputs to standardized VCF files.
- **[clamms_bed_to_vcf.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/clamms_bed_to_vcf.py)**: Converts CLAMMS BED/TSV outputs to VCF files.
- **[vcf_to_bed.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/vcf_to_bed.py)**: Extract genomic intervals from VCF to BED format.
- [vcf_filter_merge.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/vcf_filter_merge.py): Merges callsets and applies quality filters.
- [normalise_cnv_caller_quality_scores.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/normalise_cnv_caller_quality_scores.py): Normalizes caller-specific raw confidence scores to a standardized range $[0, 1]$.
- [stratify-vcf-by-acmg.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/stratify-vcf-by-acmg.py): Stratifies cohort VCF variants across ACMG class thresholds (1-5, NA, 4-5 pathogenic).

### 3. Machine Learning (ML), Feature Extraction & Quality Analytics
- [feature_extraction.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/feature_extraction.py): Ingests merged CNV callsets and extracts genomic, overlap, cohort frequency, and caller-specific features to construct ML datasets.
- [build_cohort_feature_database.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/build_cohort_feature_database.py): Compiles SQLite/DuckDB cohort-wide feature databases.
- [train_predict_catboost.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/train_predict_catboost.py): In silico prediction and filtering of candidate VCFs using trained CatBoost models with native bgzip and tabix indexing.
- [assess_callset_quality_metrics.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/assess_callset_quality_metrics.py): Multi-stage callset quality metrics extractor with 3-way comparative matrix (Total vs DEL vs DUP).
- [evaluate_caller_performance.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/evaluate_caller_performance.py): Evaluates calling accuracy against a gold standard dataset ($F_2$, MCC, Sensitivity, Precision).
- [compare_pipeline_benchmarks.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/compare_pipeline_benchmarks.py): Statistical significance and McNemar/Cochran's Q testing for pipeline comparisons.
*(Note: Offline model exploration and benchmarking scripts such as `train_xgboost.py`, `train_compare_models.py`, `evaluate_random_forest_performance.py`, and `evaluate_all_combinations.py` are located in [`benchmarks/model_comparisons/`](file:///home/ywolberg/CNV-Caller-Combiner/benchmarks/model_comparisons/).)*

### 4. Production Plotting & Performance Visualization
- **[plot_cnv_qc_metrics.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/plot_cnv_qc_metrics.py)**: Visualizes general coverage and sequencing QC metrics.
- **[plot_combination_performance.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/plot_combination_performance.py)**: Plots Venn diagrams and overlapping metrics of multi-caller combinations.
- **[plot_evaluation_performance.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/plot_evaluation_performance.py)** & **[plot_evaluation_performance_dynamic.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/plot_evaluation_performance_dynamic.py)**: Generate precision-recall and sensitivity ROC curves for callers.
- **[plot_caller_upset.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/plot_caller_upset.py)**: UpSet plots for caller concordance and overlap sets.
- **[render_cnv_infographics.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/render_cnv_infographics.py)**: High-resolution publication infographics and Pareto frontier diagrams.
*(Note: Ad-hoc and historical plotting utilities are maintained in [`benchmarks/visualization/`](file:///home/ywolberg/CNV-Caller-Combiner/benchmarks/visualization/).)*

### 5. Containers & Environment Setup
- Docker and Singularity definitions are organized in [`containers/`](file:///home/ywolberg/CNV-Caller-Combiner/containers/) with backward-compatibility symlinks preserved in `bin/`.
