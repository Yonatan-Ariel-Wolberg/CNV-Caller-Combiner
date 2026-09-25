# Machine Learning Training & Prediction (CatBoost & XGBoost)

## 1. Methodology & Algorithm
Gradient-boosted decision trees represent one of the most effective non-linear classification approaches for genomic variant quality classification. **CNV-Caller-Combiner** utilizes **CatBoost** (by default) and **XGBoost** to distinguish genuine germline copy-number events from false-positive sequencing artifacts:
* **Overcoming Single-Caller Biases**: While individual callers rely on specific distributional assumptions (e.g. Poisson, negative binomial, or beta-binomial) that often fail in regions with extreme GC content or segmental duplications, the gradient booster learns complex multivariate boundaries between true variants and artifacts.
* **Handling Categorical & Numerical Features**: CatBoost natively handles categorical combinations (caller subsets) without target leakage, using ordered boosting to mitigate overfitting.
* **Stratified Cross-Validation**: Implements k-fold cross-validation stratified by variant size class and copy number type (`DEL` vs `DUP`) to prevent data leakage and ensure generalized performance.
* **Probabilistic Scoring**: Computes a continuous prediction probability score ($0.0 \le P \le 1.0$) representing confidence that a call is a verified true positive. Variants with $P \ge 0.50$ (or custom threshold) are prioritized.

---

## 2. Nextflow Implementation & Processes
The ML training and prediction modules are encapsulated in `modules/ml/modules-train.nf` and `modules/ml/modules-predict.nf`:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `TRAIN_CATBOOST_MODEL` | `catboost_train` | Fits the gradient booster on cohort feature matrices joined with gold-standard truth labels (`--truth_labels`). |
| `EVALUATE_MODEL_PERFORMANCE` | `catboost_train` | Computes precision, recall, F1-score, PR-AUC, and ROC curves across validation folds. |
| `RUN_CATBOOST_PREDICT` | `catboost_train` | Applies the trained model to unlabelled feature matrices, outputting ranked candidate CNVs with prediction scores. |
| `FILTER_BY_ML_SCORE` | `python` | Annotates and filters consensus VCFs using the ML confidence score. |

---

## 3. Truth-Label Input Specifications
When training a new model via `--workflow train`, a tab-delimited truth table must be supplied:
```tsv
sample_id	chrom	start	end	cnv_type	truth_label
SAMPLE_001	chr1	1000000	1050000	DEL	1
SAMPLE_001	chr2	2000000	2050000	DUP	0
```
* `truth_label = 1`: Verified true positive (e.g., confirmed via qPCR, microarray, or high-confidence orthogonal calling).
* `truth_label = 0`: Verified false positive / artifact.

---

## 4. Key Parameters & Thresholds
* `params.ml_classifier`: Classifier architecture (`'catboost'` or `'xgboost'`).
* `params.ml_score_threshold`: Classification probability cutoff (default: `0.50`).
* `params.ml_pathway`: Order of operations in full workflow:
  - `'A'`: Filter dataset-internal AF $\le 5\%$ first, then run ML prediction.
  - `'B'`: Run ML prediction first, then filter dataset-internal AF $\le 5\%$.

---

## 5. Output Files
* `ml_models/catboost_cnv_model.cbm`: Serialized production model file.
* `ml_models/model_feature_importance.tsv`: Ranked ranking of predictive features.
* `ml_predictions/<SAMPLE_ID>_ml_scored.vcf.gz`: Consensus VCF annotated with `ML_SCORE` and `ML_PASS` tags.
