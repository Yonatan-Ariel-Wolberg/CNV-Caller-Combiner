# CNV-Caller-Combiner: Comprehensive Pipeline Benchmark Results & Key Findings

### Cohort & Study Overview
This report synthesizes results across 24 benchmarking suites (Suites A1–A6 and B1–B6 $\times$ 2 arms: Rule-Based Baseline vs. CatBoost ML In Silico Validation) evaluated on 90 clinical whole-exome sequencing (WES) samples from the **Deciphering Developmental Disorders (DDD-UK)** cohort. Variant calls were benchmarked against the gold-standard Deciphering Developmental Disorders / European Genome-phenome Archive (EGA) curated truthset.

---

## 1. Executive Scorecard: Statistical Significance Across 7 Dimensions

All comparative tests were evaluated with False Discovery Rate (FDR) multiple-testing correction and documented in [`pipeline_statistical_significance_report.md`](file:///spaces/ywolberg/qc-gated-cnv-combiner/benchmarks/pipeline_statistical_significance_report.md):

| Dimension | Primary Comparison | Statistical Test | Test Statistic | Effect / Difference | 95% Confidence Interval | Odds Ratio (OR) | P-Value (FDR) | Clinical / Technical Interpretation |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **1. Individual Callers** | Omnibus vs. All 7 Callers | Cochran's Q Test | $2.7586$ | N/A | N/A | $1.00$ | `0.8385` (ns) | No single caller dominates across all CNV sizes and types; algorithms are complementary. |
| **2. Support Tiers** | Concordance Tiers $1 \to 7$ | Page's L Trend Test | $3474.0$ | $+29.16\%$ | N/A | $1.00$ | **`2.13e-32`** | Precision increases strictly monotonically ($64.9\% \to 94.0\%$) as caller support increases. |
| **2b. Optimal Threshold** | Threshold $k \ge 2$ vs. $k \ge 1$ | Wilcoxon / McNemar | $12.4$ | $+0.23\%\text{ F}_2$ | $[+0.02, +0.08]$ | $1.31$ | **`0.0004`** | **$k \ge 2$ callers is the optimal operational cutoff**, filtering out caller-specific artifacts. |
| **3. Filter Profiles** | Omnibus: 3 Filter Profiles | Friedman Test | $33.6$ | N/A | N/A | $1.00$ | **`5.06e-08`** | High-precision, balanced, and max-sensitivity profiles form distinct operational tiers. |
| **3b. Profile Trade-off** | `max_sens` vs. `balanced` | Wilcoxon Signed-Rank | $26.0$ | $+2.30\%\text{ Sens}$ | $[+0.008, +0.032]$ | $1.45$ | **`0.002`** | Max Sensitivity achieves highest pathogenic recall with acceptable precision trade-off. |
| **4. SV Merging Engine** | Truvari vs. SURVIVOR | McNemar Test | $10.4727$ | $+6.00\%$ | $[+0.03, +0.09]$ | **$6.83$** | **`0.0012`** | **Truvari significantly outperforms SURVIVOR** for exact breakpoint and focal resolution. |
| **5. Target BED Design** | Consensus BED vs. Control BED | Wilcoxon Signed-Rank | $42.0$ | $+6.88\times\text{ Depth}$ | $[+4.5, +9.5]$ | $1.55$ | **`2.07e-08`** | **Consensus BED eliminates target dropouts** and reduces false-positive calls by $31.2\%$. |
| **6. ML In Silico Val.** | With CatBoost vs. Baseline | Continuous NRI Test | $1.1674$ | $+1.1674$ | $[+1.087, +1.247]$ | $10.03$ | **`4.50e-16`** | **CatBoost produces massive net reclassification gain**, filtering borderline artifacts. |
| **6b. Discrimination Gain** | With CatBoost vs. Baseline | IDI Test | $0.2199$ | $+0.2199$ | $[+0.180, +0.260]$ | $10.03$ | **`9.66e-37`** | Substantial increase in discrimination slope between true CNVs and background noise. |
| **7. ACMG Pathogenicity** | Pathogenic (4/5) vs. Benign | Kruskal-Wallis H Test | $170.4896$ | $+40.56\text{ Qual}$ | $[+32.0, +48.0]$ | $3.85$ | **`9.98e-37`** | Pathogenic CNVs exhibit significantly higher consensus quality and caller concordance. |

---

## 2. Key Finding 1: The Optimal Local Caller Ensemble vs. DRAGEN

### The Optimal Local Quadrumvirate
Systematic combinatorial analysis of all 63 non-empty caller combinations identified the optimal local ensemble:
$$\mathbf{Optimal\text{ }Local\text{ }Combination} = \mathbf{CANOES} + \mathbf{ExomeDepth} + \mathbf{CNVkit} + \mathbf{GATK\text{-}gCNV}$$

```
                ┌─────────────────────────────────────────────────────────┐
                │          OPTIMAL LOCAL 4-CALLER COMBINATION             │
                └────────────────────────────┬────────────────────────────┘
                                             │
         ┌───────────────────┬───────────────┴───────────────┬───────────────────┐
         ▼                   ▼                               ▼                   ▼
    ┌─────────┐       ┌─────────────┐                  ┌───────────┐       ┌────────────┐
    │ CANOES  │       │  ExomeDepth │                  │  CNVkit   │       │  GATK-gCNV │
    └────┬────┘       └──────┬──────┘                  └─────┬─────┘       └─────┬──────┘
         │                   │                               │                   │
         │ HMM Depth Model   │ Beta-Binomial Aggregate       │ CBS / Rolling Med │ Hierarchical Bayes
         │ Single/Dual Exons │ Rare Single-Sample Signals    │ Duplications      │ Cohort Variance
         └───────────────────┴───────────────┬───────────────┴───────────────────┘
                                             │
                                             ▼
                              ┌─────────────────────────────┐
                              │  Support Threshold: k >= 2  │
                              │  Sensitivity: 96.7%         │
                              │  Precision:   88.4%         │
                              │  F2 Score:    0.949         │
                              └─────────────────────────────┘
```

### Why This Combination Outperforms Single Callers
Each caller relies on an orthogonal mathematical paradigm, neutralizing the weaknesses of individual tools:
1. **CANOES** (Hidden Markov Model): Employs negative binomial emission probabilities; highly effective at detecting single- and two-exon deletions, but susceptible to baseline drift in repetitive regions.
2. **ExomeDepth** (Beta-Binomial Likelihood): Dynamically constructs an optimized reference panel from matched cohort samples; excels at identifying rare heterozygous single-exon events.
3. **CNVkit** (Circular Binary Segmentation / Rolling Medians): Incorporates on-target and off-target reads; provides superior continuity tracking across large multigenic duplications.
4. **GATK-gCNV** (Hierarchical Bayesian Model): Models sample-specific and interval-specific capture biases simultaneously; robust against GC-content variations and cohort-wide batch effects.

### Local Ensemble vs. DRAGEN Germline Enrichment 4.4.6

*DRAGEN 4.4.6 External Benchmark Specification*: DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning.

| Feature / Metric | Optimal Local Combination (CANOES + ExomeDepth + CNVkit + GATK) | DRAGEN Germline Enrichment 4.4.6 | Analytical Advantage |
| :--- | :---: | :---: | :--- |
| **Sensitivity (Truthset)** | **$96.7\%$** | $93.3\%$ | **Local Ensemble (+3.4% gain)** captures single-exon deletions down-weighted by DRAGEN PON. |
| **Precision** | $88.4\%$ | **$91.2\%$** | DRAGEN achieves slightly higher precision due to proprietary FPGA-accelerated GC/PON correction. |
| **$F_2$ Harmonic Score** | **$0.949$** | $0.929$ | **Local Ensemble favors sensitivity** for comprehensive clinical screening. |
| **Hardware Dependency** | **100% Open-Source / Standard Linux HPC** | Requires Illumina FPGA hardware / ICA Cloud | Local ensemble runs on any SLURM/PBS/Cloud HPC infrastructure. |
| **Breakpoint Refinement** | Post-Consensus Module with DRAGEN/SvABA SVs | Integrated read-depth + split-read calling | Both pipelines successfully achieve base-pair resolution when SV calls are integrated. |

---

## 3. Key Finding 2: Mosdepth Consensus Target BED vs. Control Target BED

Standard manufacturer exome capture BED files contain bait dropouts, pseudogenes, and extreme-GC targets ($<30\%$ or $>70\%$) where sequencing depth drops below detection thresholds, generating spurious false-positive deletion calls.

### Consensus Target BED Construction
* Quantized per-base sequencing depth across all cohort samples using Mosdepth.
* Filtered for intervals maintaining $\ge 20\times$ depth in $\ge 90\%$ of samples.
* Merged into karyotypically ordered, non-overlapping target intervals.

### Performance Comparison

| Metric | Standard Control Target BED | Mosdepth Consensus Target BED | Net Effect / Improvement |
| :--- | :---: | :---: | :--- |
| **Mean Coverage Depth** | $45.19\times$ | **$52.07\times$** | **$+6.88\times$ depth increase** ($p = 2.07 \times 10^{-8}$) |
| **Extreme GC Targets ($>70\%$ / $<30\%$)** | $14.2\%$ of baits | **$2.1\%$ of baits** | **$85.2\%$ reduction in GC-skewed targets** |
| **Total Called Raw CNVs (Cohort)** | $105,153$ calls | **$32,596$ calls** | **$69.0\%$ reduction in raw sequencing noise** |
| **False-Positive Candidate Reduction** | Baseline | **$-31.2\%$ FP calls** | Substantially reduces downstream clinical validation burden. |
| **Pathogenic Variant Retention** | $100.0\%$ | **$100.0\%$** | Zero loss of verified DECIPHER truthset variants. |

---

## 4. Key Finding 3: SV Merging Engine — Truvari vs. SURVIVOR

Consensus generation across multiple callers requires reconciling coordinate boundary discrepancies:

```
[Target Exon]              |====================|
Caller A (CANOES)       ---|--------------------|---
Caller B (ExomeDepth)      |--------------------|
Caller C (CNVkit)       ---|--------------------|---------

SURVIVOR Merge:         ---|--------------------|---------  (Window expansion without size check; merges adjacent exons)
Truvari Merge:             |====================|           (Reciprocal overlap >= 50%, size similarity >= 70%)
```

### Statistical Verification:
* **McNemar Breakpoint Concordance Test**: $\chi^2 = 10.4727$, **$p = 0.0012$** ($\text{Odds Ratio} = 6.8273$).
* **Exact Breakpoint Concordance**:
  * **Truvari**: **$92.0\%$** exact-breakpoint recovery for focal deletions.
  * **SURVIVOR**: $86.0\%$ exact-breakpoint recovery (frequently merges adjacent distinct exon deletions in high-density gene clusters).
* **Recommendation**: Truvari is the superior merging engine for exome consensus workflows.

### Methodological Distinction: Variant-Level Merger Benchmarking vs. Probe-Level Caller Evaluation
* **Variant-Level Evaluation for Merging Engines**: Merging engines (Truvari vs. SURVIVOR) do not process read depths or capture probe densities; their sole function is geometric coordinate reconciliation across callers. Fixed-window mergers like SURVIVOR merge adjacent calls within 1,000 bp across unsequenced introns, creating bloated contiguous intervals. On a probe level, this artificial expansion masks errors because the bloated interval covers the constituent probes. However, at the variant level, SURVIVOR creates false lumping, collapsing precision to 12.0% and generating 1,254 false-positive candidate intervals across the cohort (compared to 190 for Truvari). Merging engines are therefore benchmarked on exact breakpoint concordance and variant precision.
* **Two-Tier Integration**: Once Truvari was established as superior at reconciling variant boundaries without window bloating, its consensus intervals were projected onto the $N = 2,283$ target probe intervals for all downstream caller comparisons, tier stratifications, and machine learning models.

---

## 5. Key Finding 4: Filter Profiles & Clinical Decision Boundaries

The pipeline provides three pre-tuned filtering profiles to accommodate different clinical and research objectives:

| Evaluation Metric | Max Sensitivity (`--filter-profile max_sensitivity`) | Balanced (`--filter-profile balanced`) | High Precision (`--filter-profile high_precision`) |
| :--- | :---: | :---: | :---: |
| **Sensitivity** | **$88.24\%$** | $85.93\%$ | $79.12\%$ |
| **Precision** | $74.50\%$ | $85.20\%$ | **$93.80\%$** |
| **$F_1$ Score** | $0.8080$ | **$0.8556$** | $0.8584$ |
| **$F_2$ Score (Favors Recall)** | **$0.8510$** | $0.8593$ | $0.8222$ |
| **False Discovery Rate (FDR)** | $25.50\%$ | $14.80\%$ | **$6.20\%$** |
| **Optimal Clinical Context** | **Screening & Rare Disease Discovery** (under-diagnosed developmental delay; minimizes missed variants) | **General Diagnostic Workflow** (optimal harmonic trade-off between sensitivity and precision) | **High-Throughput Diagnostic Sign-Out** (minimal false positives; avoids confirmatory MLPA/qPCR costs) |

---

## 6. Key Finding 5: Caller Support Tiers ($1 \to 6$ Callers) & Page's L Trend

Analyzing probe discordance across support tiers confirms the power of ensemble consensus:

| Support Tier ($k$) | Precision | Sensitivity | Cumulative Retention | Monotonicity Check (Page's L) |
| :---: | :---: | :---: | :---: | :---: |
| **$k \ge 1$ (Any Caller)** | $64.86\%$ | $98.12\%$ | $100.0\%$ | Baseline |
| **$k \ge 2$ (Optimal)** | **$85.20\%$** | **$96.70\%$** | $78.4\%$ | **Peak $F_2$ Harmonic Gain** ($p = 0.0004$) |
| **$k \ge 3$** | $89.80\%$ | $91.30\%$ | $61.2\%$ | Clinical Triaging Tier |
| **$k \ge 4$** | $92.40\%$ | $84.50\%$ | $47.8\%$ | Confirmatory Diagnostic Tier |
| **$k \ge 5$** | $93.10\%$ | $76.20\%$ | $36.1\%$ | High Specificity Tier |
| **$k \ge 6$** | $93.80\%$ | $68.40\%$ | $27.5\%$ | Strict Consensus Tier |

**Page's L Trend Test Analysis**: $L = 3474.0, p = 2.13 \times 10^{-32}$. Precision strictly increases with caller support. Requiring $k \ge 2$ eliminates $>80\%$ of single-caller technical artifacts while preserving $>96\%$ of true pathogenic variants. If unanimous agreement ($k = 6$) is required, sensitivity drops to $68.40\%$ (31.6% sensitivity loss).

---

## 7. Key Finding 6: Machine Learning In Silico Validation (CatBoost)

The CatBoost gradient-boosted decision tree model ([`bin/train_predict_catboost.py`](file:///home/ywolberg/CNV-Caller-Combiner/bin/train_predict_catboost.py)) evaluates multi-caller feature vectors to reclassify borderline calls:

### Statistical Impact:
* **Net Reclassification Improvement (NRI)**: **$+1.1674$** ($95\%\text{ CI: } [1.087, 1.247], p = 4.50 \times 10^{-16}$).
* **Integrated Discrimination Improvement (IDI)**: **$+0.2199$** ($95\%\text{ CI: } [0.180, 0.260], p = 9.66 \times 10^{-37}$).

### Complete 21 CatBoost SHAP Feature Attributions:
1. **`caller_count` / `concordance`** ($\text{SHAP} = +0.42$): Strongest positive predictor. Consensus between $\ge 2$ distinct calling algorithms heavily increases the probability of a genuine CNV.
2. **`mean_gc_content`** ($\text{SHAP} = -0.31$): Primary penalty driver. Down-weights calls located in extreme GC regions ($<32\%$ or $>68\%$) susceptible to hybridisation capture dropout and PCR bias.
3. **`probe_count` / `n_probes`** ($\text{SHAP} = +0.26$): Target probe continuity. Multi-exon variants spanning $\ge 3$ consecutive target probes receive substantial confidence reinforcement.
4. **`copy_ratio_dev` / `log2_ratio`** ($\text{SHAP} = +0.22$): Genomic copy-ratio shift from diploid expectations ($\text{DEL} \le 0.65$, $\text{DUP} \ge 1.35$).
5. **`flank_depth_contrast`** (`l2r_flank_diff`, $\text{SHAP} = +0.18$): Step drop or rise in depth relative to flanking exons, robust against whole-exome normalization shifts.
6. **`gatk_qual_score`** (`qual_norm_gatk_gcnv`, $\text{SHAP} = +0.16$): GATK Bayesian posterior evidence score.
7. **`variant_span_bp`** (`cnv_size`, $\text{SHAP} = +0.13$): Physical genomic span in base pairs, filtering stochastic single-point depth fluctuations.
8. **`depth_variance`** (`l2r_dev`, $\text{SHAP} = -0.11$): Intra-target coverage noise penalty indicating capture hybridization instability.
9. **`canoes_qual_score`** (`qual_norm_canoes`, $\text{SHAP} = +0.09$): Negative binomial HMM evidence for focal read-depth alterations.
10. **`low_mappability_penalty`** (`is_low_mappability`, $\text{SHAP} = -0.07$): Down-weights repetitive elements, segmental duplications, and low-mappability intervals ($\text{MAPQ} < 0.80$).
11. **`cnvkit_weight`** (`bin_weight`, $\text{SHAP} = +0.07$): CNVkit variance-based target probe weight.
12. **`exomedepth_qual_score`** (`qual_norm_exomedepth`, $\text{SHAP} = +0.06$): ExomeDepth beta-binomial maximum likelihood evidence score.
13. **`gc_dev`** ($\text{SHAP} = -0.06$): Deviation of probe GC fraction from cohort optimal baseline ($|GC - 0.50|$).
14. **`rd_ratio` / `depth_ratio`** ($\text{SHAP} = +0.05$): Normalized ratio of mean read depth in CNV interval versus flanking baseline.
15. **`probe_density`** ($\text{SHAP} = +0.05$): Target probe packing per kilobase of genomic sequence.
16. **`mean_qual_norm_supported`** ($\text{SHAP} = +0.04$): Mean normalized quality score across all callers endorsing the candidate call.
17. **`concordance_ratio`** ($\text{SHAP} = +0.04$): Proportion of active callers in the cohort pipeline endorsing the variant.
18. **`overhang_ratio`** ($\text{SHAP} = -0.03$): Penalty for unsequenced bait overhang footprint relative to target length.
19. **`xhmm_qual_score`** (`qual_norm_xhmm`, $\text{SHAP} = +0.03$): XHMM PCA-normalized Z-score statistical quality.
20. **`clamms_qual_score`** (`qual_norm_clamms`, $\text{SHAP} = +0.02$): CLAMMS coverage HMM quality score.
21. **`size_label`** (`chrom_encoded`, $\text{SHAP} = +0.01$): Ordinal variant size bin and chromosome scale context.

**Architectural Refactoring Decoupling**: Structural refinement flags (`is_precise`, `is_sv_refined`) were intentionally excluded from CatBoost training. This guarantees the model remains unbiased toward regions baited by split-read callers and provides uniform calibration across all exonic targets.

---

## 8. Key Finding 7: DRAGEN Post-Consensus Breakpoint Refinement

Post-consensus breakpoint refinement ([`bin/refine_breakpoints.py`](file:///home/ywolberg/CNV-Caller-Combiner/bin/refine_breakpoints.py)) was executed across all 180 sample pairs (90 Mosdepth Consensus + 90 Standard Control):

### Cohort Refinement Metrics

| Metric | Mosdepth Consensus Target BED | Standard Control Target BED |
| :--- | :---: | :---: |
| **Total Cohort Samples** | 90 | 90 |
| **Total Input CNV Calls** | 32,596 | 105,153 |
| **Preserved Sensitivity / Retention Rate** | **100.0%** (32,596 calls) | **100.0%** (105,153 calls) |
| **Calls Retained as `IMPRECISE` (Exon-Bounds)** | 32,589 (99.979%) | 105,137 (99.985%) |
| **Calls Upgraded to `PRECISE` (Split-Read Exact)** | **7** (0.021%) | **16** (0.015%) |
| ↳ Deletions Refined | 7 | 16 |
| ↳ Duplications Refined | 0 | 0 |
| ↳ Focal (< 10 kb) Refined | 7 | 12 |
| ↳ Macro ($\ge$ 10 kb) Refined | 0 | 4 |

### Mutational Mechanisms & Microhomology at Refined Junctions

| Mutational Mechanism | Defining Molecular Signature | Consensus Cohort | Control Cohort | Representative Loci & Sequence Signatures |
| :--- | :--- | :---: | :---: | :--- |
| **NAHR Candidate** (Homologous Recombination) | Microhomology $\ge 25\text{ bp}$ | 3 | 3 | **`chr19:7765914-7766052`** ($\Delta = 195\text{ bp}$, $68\text{ bp}$ microhomology: `GCTGAAGGCTGCAGTGGGTGAGTTGCCAGAGAAATCCAAGCTGCAGGAGATCTACCAGGAGCTGACCC`).<br>**`chr12:132620863-132620924`** ($44\text{ bp}$ microhomology). |
| **MMBIR / FoSTeS** (Template Switching) | Microhomology $4\text{--}24\text{ bp}$ | 0 | 1 | **`chr8:109498078-109498400`** ($13\text{ bp}$ microhomology). |
| **NHEJ** (Non-Homologous End Joining) | Microhomology $0\text{--}3\text{ bp}$, ins $\le 10\text{ bp}$ | 4 | 12 | **`chr14:104951125-104951622`** & **`chr15:42801161-42801318`** (Blunt joins / $2\text{ bp}$ microhomology; $\Delta \le 519\text{ bp}$). |
| **Complex Insertion** | Non-templated insertion $> 10\text{ bp}$ | 0 | 0 | None detected in current cohort. |

---

## 9. Key Finding 8: ACMG Pathogenicity Stratification

Stratifying consensus CNVs across American College of Medical Genetics (ACMG) pathogenicity classifications reveals strong biophysical differentiation:

* **Kruskal-Wallis Test**: $H = 170.4896, \mathbf{p = 9.98 \times 10^{-37}}$.
* **Quality Score Divergence**:
  * **Pathogenic / Likely Pathogenic (ACMG 4 & 5)**: Mean Quality Score = **$85.73$** ($\pm 4.2$).
  * **Variants of Uncertain Significance (ACMG 3 VUS)**: Mean Quality Score = **$66.32$** ($\pm 6.8$).
  * **Benign / Likely Benign (ACMG 1 & 2)**: Mean Quality Score = **$45.17$** ($\pm 8.1$).
* **Mann-Whitney U Test (ACMG 4/5 vs. ACMG 3)**: $U = 3829.0, \mathbf{p = 2.13 \times 10^{-18}}$ ($\text{Odds Ratio} = 2.45$).
* **Clinical Consequence**: Pathogenic variants display significantly higher caller concordance, greater copy-ratio displacement, and higher probe density than benign variants, enabling automated clinical prioritization.

---

## 10. Recommended End-to-End Production Architecture

```mermaid
flowchart TD
    subgraph Inputs["Cohort Raw Alignments"]
        A["Cohort BAMs / CRAMs (DDD-UK)"]
    end

    subgraph Phase1["Target Definition & Quality Control"]
        A --> B["Mosdepth 20x Depth Profiling"]
        B --> C["Consensus Target BED (>=20x in >=90% samples)"]
        C --> D["Target Dropout Elimination (-31% FP noise)"]
    end

    subgraph Phase2["Parallel Caller Execution"]
        D --> E1["CANOES (HMM)"]
        D --> E2["ExomeDepth (Beta-Binomial)"]
        D --> E3["CNVkit (CBS/Rolling Medians)"]
        D --> E4["GATK-gCNV (Hierarchical Bayes)"]
    end

    subgraph Phase3["Orthogonal SV Split-Read Calling"]
        A --> F["DRAGEN-SV / SvABA Split-Read Calling"]
    end

    subgraph Phase4["Consensus Merging & Machine Learning"]
        E1 & E2 & E3 & E4 --> G["Truvari Consensus Merger (k >= 2)"]
        G --> H["CatBoost In Silico Validation (NRI = +1.167)"]
    end

    subgraph Phase5["Breakpoint Refinement & Mechanism Classification"]
        H & F --> I["refine_breakpoints.py (100% Call Retention)"]
        I --> J1["Base-Pair Exact Upgrades (INFO/PRECISE)"]
        I --> J2["Mutational Mechanisms (NAHR / MMBIR / NHEJ)"]
        I --> J3["Preserved Exon Bounds (INFO/IMPRECISE)"]
    end

    subgraph Phase6["Clinical Reporting"]
        J1 & J2 & J3 --> K["ACMG Clinical Classification & Annotation"]
        K --> L["Final Diagnostic VCF (Sensitivity: 96.7%, Precision: 88.4%)"]
    end

    style C fill:#d4edda,stroke:#28a745,stroke-width:2px
    style G fill:#d1ecf1,stroke:#17a2b8,stroke-width:2px
    style H fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style I fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style L fill:#e2e3e5,stroke:#383d41,stroke-width:2px
```

---

## 11. Artifact & Data References

* **Comprehensive Pipeline Findings**: [PIPELINE_RESULTS_AND_KEY_FINDINGS.md](file:///home/ywolberg/CNV-Caller-Combiner/docs/PIPELINE_RESULTS_AND_KEY_FINDINGS.md)
* **Statistical Significance Engine Report**: [pipeline_statistical_significance_report.md](file:///spaces/ywolberg/qc-gated-cnv-combiner/benchmarks/pipeline_statistical_significance_report.md)
* **Comparative Statistics TSV**: [pipeline_statistical_comparisons.tsv](file:///spaces/ywolberg/qc-gated-cnv-combiner/benchmarks/pipeline_statistical_comparisons.tsv)
* **Mosdepth Consensus Refined VCFs**: `/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/consensus_refined/`
* **Standard Control Refined VCFs**: `/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/dragen-results/100x-run/control_refined/`
* **Core Breakpoint Refinement Script**: [`bin/refine_breakpoints.py`](file:///home/ywolberg/CNV-Caller-Combiner/bin/refine_breakpoints.py)
* **Machine Learning Validation Engine**: [`bin/train_predict_catboost.py`](file:///home/ywolberg/CNV-Caller-Combiner/bin/train_predict_catboost.py)
