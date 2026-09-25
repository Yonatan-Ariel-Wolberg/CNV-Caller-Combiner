# Benchmarking, Model Comparisons, and Visualization

This directory contains standalone benchmark utilities, model comparison scripts, and ad-hoc visualization tools. These scripts are maintained for offline analysis, model development comparisons, and publication figure generation, cleanly decoupled from the production Nextflow pipeline processes in `bin/`.

## Directory Structure

### `model_comparisons/`
- `train_compare_models.py`: Evaluates and compares different classifier backends (CatBoost, XGBoost, Random Forest, LightGBM).
- `evaluate_random_forest_performance.py`: Standalone evaluation script for Random Forest baseline models.
- `evaluate_all_combinations.py`: Exhaustive combinatoric evaluation of multi-caller unions and intersections.
- `run_significance_tests.py`: Statistical significance test driver for caller comparisons.
- `train_xgboost.py`: Standalone XGBoost training and exploration script with cross-validation, SMOTE, and SHAP interpretability.

### `visualization/`
- `generate_breakdowns.py`: Generates breakdown statistics across caller callsets.
- `generate_final_plots.py`: Generates consolidated summary charts.
- `generate_poster_and_slide_figures.py`: Generates publication/conference poster figures and high-resolution slides.
- `generate_reports.py`: Ad-hoc reporting script for performance metrics across samples.
- `plot_dragen_vs_local.py`: Compares DRAGEN cloud callset performance against local callers.
- `plot_phelelani_combination_performance.py`: Historical multi-caller combination performance plotter.
- `plot_stages.py` / `plot_stages_all_metrics.py` / `plot_stages_comparison.py` / `plot_stages_extended.py` / `plot_stages_truvari_qc.py`: Multi-stage attrition and benchmark visualization scripts.
- `plot_truvari_vs_survivor.py`: Compares Truvari vs SURVIVOR merge concordances.
- `plot_workflows_comparison.py`: Cross-workflow performance comparison plotter.
- `track_pipeline_stages.py`: Utility to monitor and log stages during execution.
