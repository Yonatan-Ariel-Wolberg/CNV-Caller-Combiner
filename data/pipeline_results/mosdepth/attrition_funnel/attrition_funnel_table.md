# Variant Attrition Funnel: Discovery to Clinical Diagnostic Yield

| Pipeline Stage | Candidate Variants | Step Attrition (%) | Retention vs. Raw (%) | Biological Purpose |
| :--- | :---: | :---: | :---: | :--- |
| **1. Raw Multi-Caller Discovery** | **46,791** | `0.0%` | `100.00%` | Multi-caller parallel baseline |
| **2. Quality Score Normalization Gate** | **67** | `−99.9%` | `0.14%` | Universal scale [0–1000] & quality filters |
| **3. Intra-Sample Multi-Caller Union** | **60** | `−10.4%` | `0.13%` | Union of all caller calls (retains >=1 caller) |
| **4. Breakpoint Refinement (Precise & Imprecise)** | **60** | `0.0%` | `0.13%` | Base-pair junctions via SvABA (100% calls retained) |
| **5. Inter-Sample Cohort Union** | **58** | `−3.3%` | `0.12%` | Union across cohort samples (retains >=1 sample) |
| **6. Internal Cohort AF Filtering (<1%)** | **9** | `−84.5%` | `0.02%` | Purges benign common cohort polymorphisms |
| **7. CatBoost Machine Learning (P ≥ 0.50)** | **1** | `−88.9%` | `0.00%` | In silico gradient boosted probability triage |
| **8. ACMG Class 4/5 Prioritisation** | **0** | `−100.0%` | `0.00%` | High-confidence clinical pathogenic candidates |
