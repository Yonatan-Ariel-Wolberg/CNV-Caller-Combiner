# Automated Testing Suite (`tests/`)

This directory contains automated integration and unit tests built with `pytest` to verify the components, config parameters, helper scripts, and Nextflow processes in the pipeline.

---

## 🧪 Test Categories

### 1. Caller Integration Tests
These tests check caller modules, parameter settings, and coordinate parsers:
- **[test_clamms_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_clamms_module.py)**
- **[test_cnvkit_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_cnvkit_module.py)**
- **[test_dragen_min_samples.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_dragen_min_samples.py)** / **[test_dragen_tool_annotation.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_dragen_tool_annotation.py)**
- **[test_gatk_gcnv_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_gatk_gcnv_module.py)**
- **[test_exomedepth_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_exomedepth_module.py)**
- **[test_xhmm_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_xhmm_module.py)**

### 2. General Helper & Conversion Tests
- **[test_vcf_conversion.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_vcf_conversion.py)**: Ensures individual caller output format mappings (CSV/TSV/BED) to VCF format align correctly.
- **[test_bam_path_derivation.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_bam_path_derivation.py)**: Validates CRAM/BAM file mapping and directory structure parsing.
- **[test_normalise_quality_scores.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_normalise_quality_scores.py)**: Verifies caller-specific quality score normalization values match boundaries.

### 3. ML & Pipeline Logic
- **[test_feature_extraction.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_feature_extraction.py)**: Asserts feature extraction constructs ML training datasets without coordinate or column mismatches.
- **[test_train_module.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_train_module.py)** / **[test_train_xgboost.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_train_xgboost.py)**: Verifies training models converge and save output checkpoints.
- **[test_main_workflow.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_main_workflow.py)**: Performs integration checks on whole pipelines.

### 4. Cluster Profile Tests
- **[test_wits_cpu_profile.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_wits_cpu_profile.py)** / **[test_wits_gpu_profile.py](file:///home/ywolberg/CNV-Caller-Combiner/tests/test_wits_gpu_profile.py)**: Ensure CPU/GPU process configurations match Slurm allocation guidelines.

---

## 🏃 Running Tests

To run the full suite of tests locally or within an allocated cluster node, execute:

```bash
pytest tests/
```
