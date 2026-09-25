# Variant Attrition Funnel: Discovery to Clinical Diagnostic Yield

| Pipeline Stage | Candidate Variants | Step Attrition (%) | Retention vs. Raw (%) | Biological Purpose |
| :--- | :---: | :---: | :---: | :--- |
| **1. Raw Multi-Caller Discovery** | **94,359** | `0.0%` | `100.00%` | Multi-caller parallel baseline |
| **2. Quality Score Normalization Gate** | **964** | `−99.0%` | `1.02%` | Universal scale [0–1000] & quality filters |
| **3. Intra-Sample Multi-Caller Union** | **398** | `−58.7%` | `0.42%` | Union of all caller calls (retains >=1 caller) |
| **4. Breakpoint Refinement (Precise & Imprecise)** | **398** | `0.0%` | `0.42%` | Base-pair junctions via SvABA (100% calls retained) |
| **5. Inter-Sample Cohort Union** | **106** | `−73.4%` | `0.11%` | Union across cohort samples (retains >=1 sample) |
| **6. Internal Cohort AF Filtering (<1%)** | **33** | `−68.9%` | `0.03%` | Purges benign common cohort polymorphisms |
| **7. CatBoost Machine Learning (P ≥ 0.50)** | **7** | `−78.8%` | `0.01%` | In silico gradient boosted probability triage |
| **8. ACMG Class 4/5 Prioritisation** | **0** | `−100.0%` | `0.00%` | High-confidence clinical pathogenic candidates |
