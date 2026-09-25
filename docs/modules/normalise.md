# Normalization & Quality Scoring Module

## 1. Methodology & Algorithm
Individual CNV callers report vastly divergent, incomparable confidence metrics:
- **CANOES**: Reports negative-binomial MLCNV scores (typically 0–100+).
- **CLAMMS**: Reports likelihood-ratio test statistics (log-odds).
- **XHMM**: Reports Phred-scaled `SQ` scores (0–99).
- **CNVkit**: Reports log2 coverage depth ratios and weight-sum scores.
- **GATK gCNV**: Reports Phred-scaled genotype posterior probabilities (`QS` / `QSS`).
- **ExomeDepth**: Reports Bayes Factors (`BF`, often ranging from 0 to 1,000+).
- **DRAGEN**: Reports custom Phred-scaled confidence metrics.

The **Normalization Module** (`modules/normalise/normalise-quality.nf` and `modules/normalise/modules-normalise.nf`) harmonizes these disparate scales into standardized **Phred-scaled Quality Scores (0–99)**:
* **Caller-Specific Transformations**: Applies calibration curves and log-likelihood transformations tailored to each caller's statistical distribution.
* **Standardized VCF INFO Fields**: Sets standard `QUAL` fields and adds tool-prefixed attributes (e.g. `ORIG_QUAL`, `NORM_QUAL`).
* **Filtering Profiles**: Provides 3 pre-calibrated profiles:
  - **Balanced Profile** (Default): Calibrated for optimal F1-score balance between sensitivity and precision.
  - **High-Precision Profile**: High quality score cutoffs, ideal for studies with downstream wet-lab qPCR or MLPA validation constraints.
  - **Max-Sensitivity Profile**: Permissive quality thresholds, designed for rare disease discovery to prevent false negatives.

---

## 2. Nextflow Implementation & Processes
Encapsulated in `modules/normalise/normalise-quality.nf` and `modules/filters/quality-filters.nf`:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `NORMALISE_CNV_QUALITY_SCORES` | `python` | Applies mathematical transformations to map raw caller metrics to standardized Phred scores (0–99). |
| `FILTER_*_QUALITY` | `bcftools` | Evaluates normalized scores against profile thresholds, applying `PASS` or filtering tags. |

---

## 3. Profile Threshold Reference

| Metric / Filter | Balanced Profile | High-Precision Profile | Max-Sensitivity Profile |
| :--- | :---: | :---: | :---: |
| **Minimum Phred QUAL** | `30` | `60` | `15` |
| **Minimum Supporting Probes** | `3` | `5` | `2` |
| **Minimum Event Size (bp)** | `1,000` | `2,500` | `500` |
| **ExomeDepth Min BF** | `10.0` | `30.0` | `5.0` |

---

## 4. Output Files
* `out_<CALLER>/normalised/<SAMPLE_ID>.<caller>.normalised.vcf.gz`: VCF with standardized Phred-scaled `QUAL` metrics.
* `out_<CALLER>/filtered/<SAMPLE_ID>.<caller>.filtered.vcf.gz`: Quality-gated callset ready for consensus merging.
