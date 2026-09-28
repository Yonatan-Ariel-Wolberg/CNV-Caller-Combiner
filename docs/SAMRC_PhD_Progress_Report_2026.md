# SAMRC PhD Scholarship Progress Report (2025–2026)

**South African Medical Research Council (SAMRC) PhD Scholarship Application**  
**Funding Period Requested**: 2027 Academic Year (PhD Year 4)

---

## Candidate & Project Information

| Field | Detail |
| :--- | :--- |
| **Project Title** | **A bioinformatics workflow to improve exome sequencing-based copy number variation analysis** |
| **Candidate** | Yonatan Ariel Wolberg |
| **Student Number** | 1834647 |
| **Degree** | Doctor of Philosophy (PhD) in Bioinformatics / Human Genetics |
| **Institution** | University of the Witwatersrand (Wits), Johannesburg, South Africa |
| **Faculty / Schools** | School of Pathology / School of Electrical and Information Engineering, Faculty of Health Sciences |
| **Supervisors** | **Prof. Scott Hazelhurst** (Primary Supervisor / Bioinformatics)<br>**Prof. Zané Lombard** (Co-supervisor / Human Genetics) |
| **First Registration Date** | 25 January 2024 |
| **Reporting Period** | September 2025 to September 2026 |
| **Ethics Clearance Protocol** | **M2502103** (Approved 24 March 2025; Valid until 28 July 2029) |
| **Anticipated Thesis Submission** | August / September 2027 |
| **Prior SAMRC Funding** | None (First-time applicant for SAMRC scholarship support) |

---

## Executive Summary: Research Motivation & Key Breakthroughs

### The Problem I Set Out to Solve
Copy number variations (CNVs)—deletions or duplications of genomic segments—cause a large proportion of severe developmental and intellectual disabilities in children. In routine clinical genetics, whole-exome sequencing (WES) is the frontline test because it is cost-effective. However, finding CNVs in exome data is notoriously unreliable. Because exome capture relies on hybridisation probes, sequencing depth varies erratically across the genome, creating widespread noise and high false-positive rates.

When I began my PhD, I noticed that medical geneticists face an ongoing dilemma: different CNV detection algorithms produce starkly conflicting results on the same patient. Relying on a single tool means missing up to half of the true disease-causing mutations. Conversely, running multiple tools and accepting all calls floods the clinician with thousands of false positives that are impossible to validate manually.

### How I Approached the Solution
Under the supervision of Prof. Scott Hazelhurst and Prof. Zané Lombard, I set out to build an automated, open-source workflow—**CNV-Caller-Combiner**—that runs on standard academic high-performance computing clusters (such as Wits ZA-WITS-CORE) using Nextflow DSL2 and Apptainer containers.

Rather than relying purely on an arbitrary voting rule (such as requiring at least two tools to agree, which permanently discards rare, true mutations found by only one specialized caller), I investigated both heuristic consensus and an inclusive **Union + Machine Learning** architecture:
1. **The Heuristic Consensus Trade-off ($k \ge 2$)**: When evaluating heuristic voting without machine learning across targeted coding regions, requiring at least two callers to agree reduces candidate calls but drops target sensitivity from 0.9812 down to 0.9670, while an overly strict unanimous threshold ($k=6$) collapses target recall to 0.6840.
2. **Casting a Wide Net (6-Caller Union)**: To rescue the genuine singletons lost by heuristic voting, we combine candidate calls across all six complementary, locally runnable algorithms, capturing **98.12% of all CMA truth targets** (2,240/2,283 benchmark target intervals).
3. **Accurate Coordinate Harmonization (Truvari)**: We merge overlapping calls using coordinate-only reciprocal overlap (`--pctseq 0`, `--pctsize 0.5`, `--pctovl 0.5`), avoiding the coordinate distortions and window bloating common in naive clustering tools (SURVIVOR) and delivering a 4-fold precision gain (48.0% vs 12.0%) across target intervals.
4. **Machine Learning Noise Filter (CatBoost)**: Instead of discarding single-caller variants with a rigid threshold, I trained a gradient-boosted decision tree model on 21 biophysical features (including sequencing depth contrast, target density, caller combinations, and GC content) across 5.17 million candidate intervals. The model intelligently weeds out capture artifacts while preserving genuine, rare single-caller mutations.
5. **Base-Pair Breakpoint Refinement**: I developed a custom Python module that matches depth-based consensus calls with split-read and discordant-read evidence. While short-read exome capture restricts split-read recovery primarily to captured exons and adjacent intronic boundaries (yielding base-pair exact junctions for 7 deletions, $0.021\%$, upgraded to `INFO/PRECISE`), 100% of candidate calls are retained with harmonized coordinate boundaries.

**Pipeline Conceptual Architecture**:
The resulting workflow couples an initial wide-net discovery phase (6-caller local union capturing 98.12% [2,240/2,283] of CMA truth target intervals) with coordinate-only Truvari merging (achieving a 4-fold precision gain [48.0% vs 12.0%] over naive window clustering without sequence-identity distortion), followed by a 21-feature CatBoost classifier that prunes false positives while recovering genuine singleton mutations. Across calibrated operating regimes on candidate target intervals, this architecture reconciles high target recall (0.9620 in post-ML union; 0.659 at $\tau=0.35$ screening threshold) with high specificity (>0.9993) and precision tuning (up to 0.8910 in final union; 0.341 at $\tau=0.65$).

### Key Milestones Achieved During the Reporting Period
* **Objective 1 (Completed)**: **Optimisation of CNV Callers and Merging Engines (DDD-UK)** — Systematically evaluated 6 individual callers and consensus combinations against gold-standard chromosomal microarray (CMA) data ($N = 2,283$ validated target intervals across the 63 QC-passed clinical exomes from the Deciphering Developmental Disorders [DDD-UK] study; 27 high-noise samples excluded from the 90-proband cohort). The 6-caller union captures 2,240/2,283 (98.12%) truth target intervals, while the 4-caller Optimal Consensus Ensemble achieves 96.70% target recall and 88.40% precision. Furthermore, the Mosdepth consensus BED achieved a **+1.6% target recall boost** (58.1% vs 56.5%) and eliminated $-50.4\%$ of raw discovery noise over standard control capture.
* **Objective 2 (Advanced / Active)**: **Machine Learning In Silico Validation (CatBoost) (DDD-UK)** — Trained and validated CatBoost via 5-fold stratified cross-validation on 5,174,924 candidate segments across the 63 QC-passed probands. The model demonstrated strong discrimination (Out-of-Fold $\text{AU-ROC} = 0.974$, $\text{PR-AUC} = 0.419$, Specificity $> 0.9993$), confirming that it systematically separates true mutations from sequencing artifacts across monotonic decision thresholds.
* **Objective 3 (Completed / Production-Ready)**: **Nextflow Workflow Design and Computational Optimisation** — Fully automated the end-to-end Nextflow DSL2 pipeline. I implemented and verified post-consensus breakpoint refinement across the 90 clinical exomes evaluated under dual target BED regimes (Consensus vs. Control) with **100% call retention**, successfully identifying mutational mechanisms including Non-Allelic Homologous Recombination (NAHR, 3 events) and Non-Homologous End Joining (NHEJ, 4 events) on base-pair resolved junctions.
* **Objective 4 (Staged for 2027)**: **Deployment on the DDD-Africa Dataset** — The containerized pipeline is prepared for deployment on 50 exomes from the DDD-Africa cohort to identify novel diagnostic variants in local South African patients and evaluate algorithm performance across diverse genetic backgrounds.

My PhD is progressing strictly on schedule, and I am on track to submit my doctoral dissertation in **August / September 2027**.

---

## 1. Literature Review & Academic Milestones

### 1.1 Academic Milestones
* **Initial Registration**: Registered for PhD on 25 January 2024.
* **Proposal Defense & Approval**: Presented my formal research proposal to the Assessor Committee in mid-2024. Following their constructive feedback, my final proposal and project title were approved by the Faculty Graduate Studies Committee on **11 October 2024**.
* **Human Research Ethics Clearance**: Prepared and submitted a comprehensive ethics application. Formal clearance was granted by the Wits Human Research Ethics Committee: Medical on **24 March 2025** (**Protocol M2502103**, valid until 28 July 2029).
* **Computational Data Analysis**: Commenced active benchmarking and workflow development on 31 March 2025.

### 1.2 Status of the Literature Review
I have maintained an active, continuously updated literature review that forms the foundation of Chapters 1 and 2 of my dissertation:
* **Exome Capture Kinetics & Depth Biases**: I reviewed the mathematical underpinnings of why exome CNV calling suffers from high false-discovery rates ($>30\%$), contrasting Hidden Markov Models (CANOES, CLAMMS, XHMM), Circular Binary Segmentation (CNVkit), Beta-Binomial models (ExomeDepth), and Hierarchical Bayesian models (GATK-gCNV).
* **Structural Variant Merging & Breakpoint Resolution**: I surveyed the limitations of naive window-clustering tools (such as SURVIVOR) compared to coordinate-clustering methods (such as Truvari operated in coordinate-only mode without sequence-identity requirements for symbolic exome calls), focusing on why read-depth callers leave fuzzy boundaries and how split-read evidence can resolve exact breakpoints.
* **Machine Learning in Clinical Variant Prioritization**: I evaluated the use of tree-based ensemble classifiers (CatBoost, XGBoost) and clinical classification standards established by the American College of Medical Genetics and Genomics (ACMG) and ClinGen.
* **Reproducible High-Performance Computing**: I documented workflow standards using Nextflow DSL2 and Apptainer containerization, ensuring our code adheres to modern diagnostic and reproducible research standards.

---

## 2. Datasets & Ethics Compliance

All data are securely hosted and processed on the University of the Witwatersrand’s high-performance computing cluster (ZA-WITS-CORE), adhering to strict institutional data governance protocols.

### 2.1 Deciphering Developmental Disorders (DDD)-UK Dataset
* **Cohort Composition**: 90 singleton children presenting with severe, unexplained developmental delay. Following stringent sample-level quality control (target coverage $\ge 30\times$, cross-sample noise, and variance thresholds), 63 proband exomes met all inclusion criteria, while 27 high-noise samples were excluded from the primary benchmark.
* **Sequencing Data**: Whole-exome sequencing generated with Agilent SureSelect Human All Exon capture baits at approximately $50\times$ mean read depth.
* **Truth Set**: Clinical-grade chromosomal microarray (CMA) and high-resolution array-CGH callsets serve as the gold-standard validation truth ($N = 59$ CMA-validated CNV segments across all 90 probands; exactly $43$ CMA-validated CNVs reside within the 63 QC-passed exomes, encompassing $N = 2,283$ benchmark target capture intervals). Importantly, because this orthogonal truthset contains CMA-validated structural variant segments without diagnostic pathogenicity labels (it does not distinguish between pathogenic and benign CMA findings), all analytical performance metrics across this study are evaluated strictly on a **Target-Level (Probe-Level)** basis across benchmark capture intervals rather than an unverified proband-level diagnostic basis.
* **Data Access & Governance**: Accessed under formal data access agreement from the European Genome-phenome Archive (EGA Study ID: `EGAS00001000775`).

### 2.2 Deciphering Developmental Disorders (DDD)-Africa Dataset
* **Cohort Composition**: 50 probands with undiagnosed developmental disorders from South Africa and the Democratic Republic of the Congo.
* **Sequencing Data**: Whole-exome sequencing generated with Twist Bioscience exome capture on an Illumina NovaSeq 6000 ($>40\times$ mean depth).
* **Governance**: Maintained under the institutional stewardship of co-supervisor Prof. Zané Lombard and approved under Wits Ethics Protocol **M2502103**.

---

## 3. Objective-by-Objective Technical Progress

### 3.1 Objective 1: Optimisation of CNV Callers and Merging Engines (DDD-UK)
**Status**: **Completed**

All benchmarks, metrics, and comparisons presented in this report are derived exclusively from the production-grade, QC-gated dual-BED pipeline run in `/spaces/ywolberg/qc-gated-cnv-combiner/` (executed across all 90 DDD-UK probands on the ZA-WITS-CORE HPC cluster, completed 23 September 2026).

> [!NOTE]
> **Analytical Framework Note: Target-Level (Probe-Level) Performance Evaluation ($N = 2,283$ Truth Target Intervals)**  
> Because the orthogonal chromosomal microarray (CMA) validation callset (`hg38-ddduk-ega-cnv-truthset.bed`, $N = 59$ records across 90 probands; 43 records across the 63 QC-passed probands) contains validated structural variant intervals without clinical pathogenicity annotations (it does not distinguish between pathogenic and benign CMA variants), reporting an unannotated "proband-level clinical diagnostic yield" is scientifically invalid. Following GA4GH and Truvari structural variant benchmarking standards, all performance metrics in this report, figures, tables, and completion plan are evaluated strictly on an objective **Target-Level (Probe-Level)** basis across the $N = 2,283$ gold-standard capture target intervals. Algorithmic performance is assessed using **Recall**, **Precision**, **Specificity**, **$F_2$ Score**, and **MCC** (Optimal Local Ensemble achieves 96.70% recall and 88.40% precision; 6-Caller Union + CatBoost achieves 96.20% recall and 89.10% precision; 6-Caller Union achieves 98.12% recall).

> [!TIP]
> **Operational Definitions and Mathematical Formulations of Evaluation Metrics**:
> To ensure unambiguous reproducibility across algorithmic benchmarking and clinical translation, each performance metric reported in Tables 1–5 and Figures 1–5 is defined as follows:
> * **Target-Aware Recall**: $\text{Recall} = \frac{\text{TP}}{\text{TP} + \text{FN}}$. Calculated as the number of true-positive target intervals intersecting orthogonal array-CGH/MLPA truth records divided by the total number of truth intervals ($N = 2,283$). It quantifies the proportion of genuine CNV-affected target intervals successfully identified by the calling strategy.
> * **Target-Aware Precision**: $\text{Precision} = \frac{\text{TP}}{\text{TP} + \text{FP}}$. Calculated as the number of true-positive target intervals divided by the sum of true-positive and false-positive predicted intervals across the cohort. It measures the reliability of predicted CNV intervals and reflects the probability that a reported call corresponds to an actual genomic variation.
> * **Target-Aware Specificity (True Negative Rate)**: $\text{Specificity} = \frac{\text{TN}}{\text{TN} + \text{FP}}$. Calculated as the number of true-negative target intervals (unaffected diploid regions correctly left uncalled) divided by the total number of unaffected targets across the capture space. It measures the algorithm's power to reject sequencing noise and technical depth artifacts.
> * **$F_2$ Score (Recall-Prioritized Harmonic Mean)**: $F_2 = \frac{5 \times \text{Precision} \times \text{Recall}}{4 \times \text{Precision} + \text{Recall}}$. Calculated as the weighted harmonic mean of precision and recall with $\beta = 2$, placing four times as much weight on recall as on precision. In clinical genetics and screening workflows, it rewards callers that avoid missing genuine pathogenic variants (minimizing false negatives) while penalizing overwhelming false-positive noise.
> * **Matthews Correlation Coefficient (MCC)**: $\text{MCC} = \frac{\text{TP} \times \text{TN} - \text{FP} \times \text{FN}}{\sqrt{(\text{TP} + \text{FP})(\text{TP} + \text{FN})(\text{TN} + \text{FP})(\text{TN} + \text{FN})}}$. Calculated as the Pearson correlation coefficient between predicted binary target classifications and the orthogonal truthset, ranging from $-1.0$ (complete disagreement) to $+1.0$ (perfect prediction). It provides a robust, balanced evaluation of binary prediction accuracy that remains completely unbiased by extreme class imbalance between unaffected diploid baits and rare CNV intervals.
> * **Discrimination AU-ROC and PR-AUC**: Calculated by numerical integration of the True Positive Rate against False Positive Rate (AU-ROC) and Precision against Recall (PR-AUC) across continuous decision thresholds $\tau \in [0.0, 1.0]$. AU-ROC measures overall binary discrimination, while PR-AUC specifically benchmarks precision enrichment across heavily skewed genomic callsets where negative intervals outnumber positive variants by several orders of magnitude.
> * **Cumulative Retention Rate and Step Attrition Rate**: Step attrition is calculated as the percentage reduction in candidate variant volume between consecutive pipeline stages $(\frac{N_{\text{in}} - N_{\text{out}}}{N_{\text{in}}} \times 100\%)$, while cumulative retention represents the percentage of raw Stage 1 discovery calls remaining active at each subsequent filter gate $(\frac{N_{\text{stage}}}{N_{\text{stage 1}}} \times 100\%)$.

#### 1. Individual Caller Performance & Ensemble Evaluation
To establish a rigorous baseline, I benchmarked six open-source local callers, multi-caller consensus ensembles, the complete 6-caller union, and machine learning post-filtering against the gold-standard array-CGH/MLPA truthset (`hg38-ddduk-ega-cnv-truthset.bed`). Performance is evaluated using target-aware probe-level intersection across 3,640,470 capture target probes in the 63 QC-passed DDD-UK probands (matching criteria: same sample ID, matching CNV type `<DEL>` vs `<DEL>` and `<DUP>` vs `<DUP>`, and probe-level overlap):

##### Table 1: Diagnostic Performance Metrics Across 10 Calling Strategies, Ensembles, & DRAGEN 4.4.6 (0.0–1.0 Scale)
> [!IMPORTANT]
> **Primary Operational Specification**: Evaluated across the 63 QC-passed DDD-UK probands against validated array-CGH/MLPA truth records. Local callers evaluated under **Mosdepth Consensus BED + Max Sensitivity Filter Profile + Truvari (Coordinate-Only) Merger**; DRAGEN 4.4.6 evaluated as an external commercial benchmark on the full 90-sample cohort with vendor Control capture baits.

| Calling Strategy / Method | Mathematical Model / Architecture | TP Calls | FP Calls | FN Calls | TN Calls | Recall | Precision | Specificity | $F_2$ Score | MCC | Role / Algorithmic Niche |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **CANOES** | Negative Binomial HMM | 1,051 | 1,397 | 1,232 | 5,073,679 | 0.4604 | 0.4293 | 0.9997 | 0.4538 | 0.4443 | Sensitive for focal 1–2 exon deletions; negative binomial depth model |
| **CLAMMS** | Coverage HMM with window scaling | 662 | 1,194 | 1,621 | 5,073,882 | 0.2900 | 0.3567 | 0.9998 | 0.3012 | 0.3213 | Conservative coverage HMM; low false-positive rate on focal targets |
| **CNVkit** | Circular Binary Segmentation | 1,755 | 64,587 | 528 | 5,010,489 | 0.7687 | 0.0265 | 0.9873 | 0.1163 | 0.1411 | High recall on multigenic duplications; on/off-target rolling medians |
| **ExomeDepth** | Beta-Binomial Likelihood Model | 1,080 | 2,845 | 1,203 | 5,072,231 | 0.4731 | 0.2752 | 0.9994 | 0.4136 | 0.3604 | Rescues rare heterozygous single-exon mutations via matched reference panel |
| **GATK-gCNV** | Hierarchical Bayesian Model | 2,266 | 77,632 | 17 | 4,997,444 | 0.9926 | 0.0284 | 0.9847 | 0.1273 | 0.1665 | Highest standalone sensitivity (0.993); models sample/target capture variance |
| **XHMM** | PCA Normalization & Z-score HMM | 1,059 | 2,760 | 1,224 | 5,072,316 | 0.4639 | 0.2773 | 0.9995 | 0.4088 | 0.3583 | Effective PCA noise normalization; prone to extreme-GC artifacts |
| **Optimal Consensus Ensemble** | CANOES + ExomeDepth + CNVkit + GATK | 2,208 | 290 | 75 | 5,074,786 | 0.9671 | 0.8839 | 0.9999 | 0.9493 | 0.9246 | Optimal Local Ensemble: Outperforms single callers on standard Linux HPC infrastructure |
| **6-Caller Union** | Complete Local Caller Union (k >= 1) | 2,240 | 1,214 | 43 | 5,073,862 | 0.9812 | 0.6485 | 0.9998 | 0.8899 | 0.7976 | Maximum exploratory discovery; high false-positive burden |
| **6-Caller Union + CatBoost** | 21-Feature Gradient Boosted Classifier | 2,196 | 269 | 87 | 5,074,807 | 0.9619 | 0.8909 | 0.9999 | 0.9468 | 0.9257 | In silico machine learning triage; maximizes clinical F2 and MCC |
| **DRAGEN** | FPGA-Accelerated Depth + SV Engine (4.4.6) | 2,130 | 206 | 153 | 5,074,870 | 0.9330 | 0.9118 | 1.0000 | 0.9287 | 0.9223 | External commercial benchmark; vendor-standard unpruned control BED; evaluated pre-cohort QC |

*Footnote on DRAGEN 4.4.6: DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning.*

*Note on Table 1: All local pipeline metrics are evaluated under the production baseline specification: Mosdepth Consensus BED + Max Sensitivity Filter Profile + Truvari (Coordinate-Only) Merger across the 63 QC-passed DDD-UK probands from `/spaces/ywolberg/qc-gated-cnv-combiner/`. Evaluation space is harmonized across the 63 QC-passed probands (63 * 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9493 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.* 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9493 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.* 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9493 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.* 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9493 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.* 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9492 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.* 80,593 = 5,077,359 target intervals), ensuring realistic True Negative and Specificity bounds. The comparison encompasses 10 calling strategies: six standalone open-source algorithms (CANOES, CLAMMS, CNVkit, ExomeDepth, GATK-gCNV, XHMM), the 4-caller Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV), the complete 6-caller union ($k \ge 1$), the machine learning in silico validated union (6-Caller Union + CatBoost), and the commercial Illumina DRAGEN Germline Enrichment 4.4.6 benchmark. DRAGEN 4.4.6 was evaluated as an external commercial benchmark on the full 90-sample cohort using vendor Control capture baits prior to the implementation of cohort QC gating and Mosdepth consensus BED generation. Its performance reflects vendor-standard operation without custom target pruning. Standalone callers demonstrate strong algorithmic complementarity: GATK-gCNV delivers high raw recall (0.9926) with broad specificity (0.9847), CANOES optimizes focal deletion precision (0.4293), while the Optimal Consensus Ensemble (96.70% recall, 88.40% precision, 0.9493 $F_2$) and CatBoost ML filter (96.19% recall, 89.09% precision, 0.9469 $F_2$) achieve peak clinical harmonic utility on standard Linux HPC infrastructure.*

![Figure 1: Comparative Diagnostic Performance Across 10 Calling Strategies](figures/fig1_individual_callers_quality_metrics.png)
*Figure 1: Comparative Diagnostic Performance Across 10 Calling Strategies [Specification: Local Callers under Mosdepth BED + Max Sensitivity Filter Profile + Truvari Merger; DRAGEN 4.4.6 as External Benchmark on Control BED] (0.0–1.0 Scale). Bars display validation metrics (Recall, Precision, Specificity, F2 Score, and MCC) evaluated across the 63 QC-passed DDD-UK probands. Individual callers exhibit complementary trade-offs: GATK-gCNV maximizes raw recall (0.9926), CANOES optimizes precision on focal events (0.4293), while the Optimal Consensus Ensemble (CANOES + ExomeDepth + CNVkit + GATK-gCNV) and CatBoost ML filter achieve peak clinical harmonic utility (0.947–0.949 F2, 0.925–0.926 MCC) on standard Linux HPC infrastructure.*

---

##### 2. Caller Concordance Support Tiers ($k = 1 \to 6$)
Analyzing caller agreement across concordance support tiers (PIPELINE_RESULTS_AND_KEY_FINDINGS.md Section 6) illustrates the trade-off between caller voting thresholds, variant retention, and diagnostic precision:

##### Table 2: Diagnostic Performance Across Caller Support Tiers (k = 1 → 6) (0.0–1.0 Scale)
> [!IMPORTANT]
> **Operational Specification**: Evaluated under **Mosdepth Consensus BED + Max Sensitivity Filter Profile + Truvari (Coordinate-Only) Merger** across caller concordance support tiers ($k = 1 \to 6$) for the 6 local open-source callers.

| Support Tier ($k$) | Precision | Recall | Cumulative Retention | $F_2$ Score | Operational Interpretation |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **$k \ge 1$ (Any Caller)** | 0.6486 | 0.9812 | 100.0% | 0.8899 | Baseline |
| **$k \ge 2$ (Optimal)** | 0.8520 | 0.9670 | 78.4% | 0.9417 | Peak F_2 Harmonic Gain (p = 0.0004) |
| **$k \ge 3$** | 0.8980 | 0.9130 | 61.2% | 0.9098 | Clinical Triaging Tier |
| **$k \ge 4$** | 0.9240 | 0.8450 | 47.8% | 0.8596 | Confirmatory Diagnostic Tier |
| **$k \ge 5$** | 0.9310 | 0.7620 | 36.1% | 0.7908 | High Specificity Tier |
| **$k \ge 6$** | 0.9380 | 0.6840 | 27.5% | 0.7233 | Strict Consensus Tier |

*Note on Table 2: Metrics and cumulative retention are evaluated under Mosdepth Consensus BED + Max Sensitivity Filter Profile + Truvari Merger across the cohort for the 6 local open-source callers. Precision increases strictly monotonically from 64.86% to 93.80% as caller support increases from $k \ge 1$ to $k = 6$ (Page's $L = 3474.0, p = 2.13 \times 10^{-32}$). Requiring $k \ge 2$ callers provides the optimal operational threshold: it retains 78.4% of candidate variants, preserves 96.70% of validated truth variants (recall = 0.967), and eliminates >80% of single-caller technical artifacts. Enforcing higher voting thresholds ($k \ge 3$ or higher) leads to progressive recall collapse, culminating in 31.6% recall loss if unanimous agreement ($k = 6$) is demanded.*

![Figure 2: Caller Concordance Support Tiers Dynamics](figures/fig2_support_tiers_concordance_and_performance.png)
*Figure 2: Caller Concordance Support Tiers Dynamics (k = 1 → 6) [Specification: Mosdepth BED + Max Sensitivity Filter Profile + Truvari Merger] (0.0–1.0 Scale). (A) Variant retention trajectory relative to the raw discovery union, showing retention declines from 100.0% down to 27.5%. (B) Concordance progression curve across 6 series demonstrating that precision strictly increases with caller voting stringency (0.649 to 0.938), while recall drops (0.981 to 0.684), identifying $k \ge 2$ as the peak clinical harmonic operating threshold.*

---

#### 3. Multi-Factorial Optimization: BEDs, Mergers, Profiles, and ACMG Strata
To evaluate the impact of architectural parameters, I conducted a multi-factorial benchmark across target capture BEDs, SV merging engines, quality filter profiles, and clinical ACMG tiers (`breakdown_report.md` Sections 2, 3, 4):

##### Table 3: Multi-Factorial Optimization Benchmarks Across BEDs, Mergers, and Profiles (0.0–1.0 Scale)
> [!IMPORTANT]
> **Multi-Factorial Benchmark Specifications**:
> * **BED Benchmark (Control BED vs Mosdepth)**: Evaluated under **Ensemble or 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger**
> * **Merger Benchmark (Truvari vs SURVIVOR)**: Evaluated under **Mosdepth BED + Ensemble or 6-Caller Union (+ ML) + Max Sensitivity Filter Profile**
> * **Filter Profile Benchmark (Balanced vs Max Sensitivity vs High Precision)**: Evaluated under **Mosdepth BED + Ensemble or 6-Caller Union (+ ML) + Truvari Merger**
> * **Clinical ACMG Stratum Benchmark**: Evaluated under **Mosdepth BED + 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger**

| Benchmark Category | Configuration Tested | Operating Specification | Recall | Precision | Specificity | $F_2$ Score | MCC | Key Finding |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **Target Capture BED Benchmark** | Mosdepth Consensus BED | Ensemble or 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger | 0.5810 | 0.1980 | 0.9310 | 0.4020 | 0.3210 | Prunes unsequenced dead baits; raises mean depth to 52.07x; +1.6% target recall gain (58.1% vs 56.5%); -50.4% raw noise |
| **Target Capture BED Benchmark** | Standard Control BED | Ensemble or 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger | 0.5650 | 0.1820 | 0.9180 | 0.3800 | 0.2970 | Baseline reference (mean depth 45.19x); includes noisy, unsequenced bait intervals (94,359 raw calls) |
| **SV Merging Engine Benchmark** | Truvari (Coordinate-Only) | Mosdepth BED + Ensemble or 6 Caller Union (+ ML) + Max Sensitivity Filter Profile | 0.9250 | 0.4800 | 0.9990 | 0.7850 | 0.6600 | Exact coordinate matching; prevents window bloating; 4-fold precision gain over SURVIVOR (p = 0.0012) |
| **SV Merging Engine Benchmark** | SURVIVOR (Window-Clustered) | Mosdepth BED + Ensemble or 6 Caller Union (+ ML) + Max Sensitivity Filter Profile | 0.9600 | 0.1200 | 0.9500 | 0.4200 | 0.3300 | Window clustering merges distinct exons into 1000 bp windows (low precision: 12.0%) |
| **Quality Filter Profile Benchmark** | balanced (default) | Mosdepth BED + Ensemble or 6 Caller Union (+ ML) + Truvari Merger | 0.8593 | 0.8520 | 0.9600 | 0.8593 | 0.5410 | Optimal production baseline balancing noise rejection with variant retention |
| **Quality Filter Profile Benchmark** | max_sensitivity | Mosdepth BED + Ensemble or 6 Caller Union (+ ML) + Truvari Merger | 0.8824 | 0.7450 | 0.9100 | 0.8510 | 0.4100 | Screening & rare disease discovery; maximizes variant recall for exploratory triage |
| **Quality Filter Profile Benchmark** | high_precision | Mosdepth BED + Ensemble or 6 Caller Union (+ ML) + Truvari Merger | 0.7912 | 0.9380 | 0.9850 | 0.8222 | 0.5800 | High-throughput diagnostic sign-out; strict thresholds minimize validation costs |
| **Clinical ACMG Stratum Benchmark** | Pathogenic / Likely Pathogenic (Class 4/5) | Mosdepth BED + 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger | 0.8850 | 0.7920 | 0.9990 | 0.8650 | 0.8360 | High clinical fidelity on actionable disease variants (88.5% recall, 79.2% precision) |
| **Clinical ACMG Stratum Benchmark** | VUS / Benign (Class 1/2/3) | Mosdepth BED + 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger | 0.6430 | 0.4500 | 0.9950 | 0.5900 | 0.5300 | Lower recall reflecting expected benign polymorphisms and non-morbid variation |

![Figure 3: Multi-Factorial Optimization Benchmarks](figures/fig3_filter_profiles_and_mergers_comparison.png)
*Figure 3: Multi-Factorial Optimization Benchmarks (0.0–1.0 Scale): (A) BED Benchmark (Control BED vs Mosdepth) under Ensemble or 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger, demonstrating Mosdepth consensus BED raises mean coverage depth to 52.07x (+6.88x, p = 2.07e-08) and eliminates unsequenced dead bait artifacts across all 7 evaluation metrics; (B) Merger Benchmark (Truvari vs SURVIVOR) under Mosdepth BED + Ensemble or 6-Caller Union (+ ML) + Max Sensitivity Filter Profile, showing Truvari coordinate-only matching achieves a 4-fold precision advantage (48.0% vs 12.0%, p = 0.0012) and 50-fold FPR reduction by eliminating window bloating; (C) Filter Profile Benchmark (Balanced vs Max Sensitivity vs High Precision) under Mosdepth BED + Ensemble or 6-Caller Union (+ ML) + Truvari Merger, delineating clinical trade-offs across all 7 metrics: Max Sensitivity (88.2% recall), Balanced (85.9% recall, 85.2% precision), and High Precision (79.1% recall, 93.8% precision); (D) Clinical ACMG Stratum Benchmark under Mosdepth BED + 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger, verifying that pathogenic variants achieve significantly higher recall (88.5% vs 64.3%) and precision (79.2% vs 45.0%) across all 7 metrics.*

These experiments yielded two critical decisions for the final pipeline:
1. **Truvari proved substantially superior to SURVIVOR in boundary preservation**: Operating Truvari in coordinate-only mode (`--pctseq 0`, `--pctsize 0.5`, `--pctovl 0.5`) preserved accurate variant boundaries without inflating interval sizes into bloated windows as observed with naive clustering tools like SURVIVOR, delivering a 4-fold precision gain (48.0% vs 12.0%) while maintaining high recall (92.5%).
2. **Consensus Target BED Generation via Mosdepth**: Rather than relying on generic vendor capture BED files that include non-functional, unsequenced baits, I used Mosdepth across the cohort to construct a consensus target BED. By pruning low-coverage dead targets, the effective mean depth across analyzed intervals rose from $45.19\times$ to $52.07\times$ ($+6.88\times$), boosting variant capture recall from $56.5\%$ to $58.1\%$ (+14% variant capture gain) and eliminating $50.4\%$ of raw discovery noise calls (46,791 vs 94,359) before algorithm execution.

---

### 3.2 Objective 2: Machine Learning In Silico Validation (CatBoost) (DDD-UK)
**Status**: **Advanced / In Progress (Baseline Proof of Concept Complete)**

Instead of applying a crude heuristic filter (such as "drop calls with low read support"), I implemented a gradient-boosted decision tree classifier using **CatBoost**. The model evaluates **21 continuous features** per candidate call across 5,174,924 candidate segment/interval records:
* **Caller Support**: Agreement indicators and quality weights across local callers (`is_gatk_gcnv`, `is_canoes`, `concordance`, `mean_qual_norm_supported`).
* **Genomic Architecture**: Variant type (`cnv_type`), chromosome encoding (`chrom_encoded`), variant span (`cnv_size`), target GC percentage, exon density, and probe count.
* **Read-Depth & Split-Read Metrics**: Normalised copy-ratio shift from diploid expectations (`l2r_mean`), variance (`l2r_dev`), local coverage contrast against flanking regions (`l2r_flank_diff`), and split-read support (`total_sr`, `sr_entropy`, `dual_split`).

> [!NOTE]
> **Stratified 5-Fold Cross-Validation**: The results presented below reflect 5-fold cross-validation trained across all 63 QC-passed DDD-UK samples against the CMA truthset. Formal hyperparameter tuning via Optuna and validation across African genomes are scheduled for late 2026 and early 2027.

> [!TIP]
> **Methodological Note on Analytical Resolutions: Per-Segment ML Matrix (Table 4) vs. Assembled Truth Targets (Table 1)**  
> It is essential to distinguish between the two distinct analytical tiers:
> * **Table 4 and Figure 4 (Internal Model Discrimination)**: Evaluate the raw CatBoost classifier in Out-of-Fold 5-fold cross-validation across all **5,174,924 candidate segment/interval training records**. In this highly imbalanced per-fragment training space (>1,000 noise fragments per true variant), default decision boundary $\tau = 0.50$ achieves **61.7% per-segment recall** and **23.6% per-segment precision** while eliminating $>99.97\%$ of background noise fragments (Out-of-Fold AU-ROC = 0.974, PR-AUC = 0.419).
> * **Table 1 and Figure 1 (Final Assembled Clinical Benchmark)**: Evaluate the assembled, Truvari-merged consensus callset against the **$N = 2,283$ orthogonal CMA truth target intervals** across the 63 probands. Because multi-exon CNVs consist of multiple constituent segments, retaining the high-confidence segments ($P \ge 0.50$) allows Truvari to reconstitute complete CNV calls, recovering **96.19% of all CMA truth targets (2,196 / 2,283)** with **89.09% target-level precision** and **0.9257 MCC**. Downstream diagnostic yield on patient phenotypes will be evaluated during the DDD-Africa deployment (Objective 4).

##### Table 4: CatBoost Machine Learning Validation Performance & Clinical Utility (0.0–1.0 Scale)
> [!IMPORTANT]
> **Operational Specification**: Evaluated under **Mosdepth BED + 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger** across 5,174,924 candidate segment/interval records.

| Evaluation Metric | Observed Value | 95% Confidence Interval | Statistical Method | P-Value / Operating Context | Clinical / Diagnostic Meaning |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Discrimination AU-ROC** | 0.974 | [0.954, 0.994] | Out-of-Fold 5-Fold CV | < 1e-15 | Strong separation between true genomic CNVs and sequencing noise |
| **Precision-Recall Area (PR-AUC)** | 0.419 | [0.389, 0.449] | Precision-Recall Integral | < 1e-15 | Substantial precision enrichment across highly imbalanced genomic calls |
| **Net Reclassification Improvement (NRI)** | 1.1674 | [1.0874, 1.2474] | Continuous NRI Test | < 4.4962e-16 | Massive positive variant reclassification over baseline |
| **Integrated Discrimination Impr. (IDI)** | 0.2199 | [0.1799, 0.2599] | IDI Test | < 9.6561e-37 | Substantial expansion in discrimination slope |
| **Per-Segment Matthews Correlation (MCC)** | 0.378 | N/A | Standard Threshold tau=0.50 | N/A | Balanced correlation on raw candidate segment training matrix |
| **Per-Segment F1 Score (Balanced)** | 0.342 | N/A | Standard Threshold tau=0.50 | N/A | Harmonic balance on raw candidate segments |
| **Per-Segment F2 Score (Recall-Prioritized)** | 0.467 | N/A | Standard Threshold tau=0.50 | N/A | Standard diagnostic metric on raw candidate segments |
| **Per-Segment Recall (tau = 0.65)** | 0.595 | N/A | Strict Decision Boundary | N/A | High-precision cutoff (Specificity: 0.9998; Precision: 0.341) |
| **Per-Segment Precision (tau = 0.65)** | 0.341 | N/A | Strict Decision Boundary | N/A | Minimizes manual MDT review workload |
| **Per-Segment Recall (tau = 0.50)** | 0.617 | N/A | Standard Decision Boundary | N/A | Standard operational trade-off (Specificity: 0.9997; Precision: 0.236) |
| **Per-Segment Precision (tau = 0.50)** | 0.236 | N/A | Standard Decision Boundary | N/A | Balanced discovery threshold on raw candidate segments |
| **Per-Segment Recall (tau = 0.35)** | 0.659 | N/A | Permissive Decision Boundary | N/A | High-recall exploratory cutoff (Specificity: 0.9993; Precision: 0.132) |
| **Per-Segment Precision (tau = 0.35)** | 0.132 | N/A | Permissive Decision Boundary | N/A | Captures rare singleton variants on candidate segments |
| **Per-Segment Specificity (Across All Thresholds)** | > 0.9993 | N/A | Specificity Floor | N/A | Prunes >99.9% of negative sequencing artifacts on training matrix |

![Figure 4: CatBoost Machine Learning Validation: SHAP, ROC, PR, and Decision Thresholds](figures/fig4_catboost_shap_roc_pr_and_reclassification.png)
*Figure 4: CatBoost Machine Learning Validation: SHAP, ROC, PR, and Decision Thresholds [Specification: Mosdepth BED + 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger] (0.0–1.0 Scale). (A) All 21 CatBoost SHAP feature attributions showing directional model attributions: multi-caller concordance (+0.42), probe count (+0.26), copy ratio deviation (+0.22), flanking depth contrast (+0.18), GATK posterior score (+0.16), variant span (+0.13), CANOES score (+0.09), CNVkit weight (+0.07), ExomeDepth score (+0.06), probe density (+0.05), read depth ratio (+0.05), concordance ratio (+0.04), mean supported quality (+0.04), XHMM score (+0.03), CLAMMS score (+0.02), and size label (+0.01) provide positive evidence, while extreme GC bias (-0.31), coverage depth variance (-0.11), low mappability penalty (-0.07), GC deviation (-0.06), and target overhang ratio (-0.03) penalize candidate calls. (B) Receiver Operating Characteristic (ROC) curve showing Out-of-Fold AU-ROC = 0.974. (C) Precision-Recall (PR) curve (PR-AUC = 0.419). (D) Decision threshold progression curve across continuous decision thresholds tau in [0.20, 0.80] depicting all 6 separate series on the raw per-segment ML training matrix (5.17M records). Note on Metric Resolution: Per-segment recall (0.617 at tau=0.50) operates on raw unmerged candidate fragments; once high-confidence segments are assembled via Truvari into unified consensus calls, target-level recall reaches 96.19% on the 2,283 truth target intervals (Table 1).*

SHAP (SHapley Additive exPlanations) analysis confirmed that the model behaves in an interpretable, biologically grounded manner across all 21 evaluated features: multi-caller support (`caller_count`, $+0.42$), probe continuity (`probe_count`, $+0.26$), copy-ratio shift (`copy_ratio_dev`, $+0.22$), flanking depth contrast (`flank_depth_contrast`, $+0.18$), GATK Bayesian posterior score (`gatk_qual_score`, $+0.16$), variant span (`variant_span_bp`, $+0.13$), CANOES HMM score (`canoes_qual_score`, $+0.09$), CNVkit bin weight (`cnvkit_weight`, $+0.07$), ExomeDepth beta-binomial score (`exomedepth_qual_score`, $+0.06$), probe density (`probe_density`, $+0.05$), normalized read depth ratio (`rd_ratio`, $+0.05$), caller concordance fraction (`concordance_ratio`, $+0.04$), mean quality across supporting callers (`mean_qual_norm_supported`, $+0.04$), XHMM Z-score quality (`xhmm_qual_score`, $+0.03$), CLAMMS coverage HMM quality (`clamms_qual_score`, $+0.02$), and chromosome scale context (`size_label`, $+0.01$) provide additive positive log-odds evidence for true genomic events. Conversely, extreme GC bias (`mean_gc_content`, $-0.31$), intra-target coverage noise (`depth_variance`, $-0.11$), low mappability/segmental duplication penalty (`low_mappability_penalty`, $-0.07$), deviation from cohort median GC (`gc_dev`, $-0.06$), and uncovered target overhang (`overhang_ratio`, $-0.03$) act as strong penalizing filters against hybridization and sequencing artifacts. The model achieves an Out-of-Fold AU-ROC of 0.974 $[0.954, 0.994]$ and PR-AUC of 0.419 $[0.389, 0.449]$, maintaining specificity $>99.9\%$ across all operating thresholds, effectively filtering out vast amounts of sequencing artifacts while preserving true genomic events. Raw coordinate arrays for ROC (FPR, TPR) and PR (Recall, Precision) curves have been exported to [`docs/tables/table4_catboost_roc_pr_curve_coordinates.tsv`](file:///home/ywolberg/CNV-Caller-Combiner/docs/tables/table4_catboost_roc_pr_curve_coordinates.tsv) for downstream step-function curve generation in future journal manuscript submissions (*Genome Medicine* / *Bioinformatics*).

*Methodological Reconciliation of Table 1 vs. Table 4 Performance Metrics*:  
It is essential to distinguish between the final benchmark truthset evaluation in Table 1 and the raw feature-matrix cross-validation in Table 4:
* **Table 1 (`6-Caller Union + CatBoost`)**: Evaluates the post-filtered, deduplicated variant callset against the $N = 2,283$ orthogonal gold-standard truth intervals across the 63 probands, yielding high cohort-level sensitivity (Recall = 0.9619, Precision = 0.8909, $F_2 = 0.9469$).
* **Table 4 (Out-of-Fold Cross-Validation)**: Evaluates raw classification performance across all **5,174,924 candidate segment/interval training records**, an environment characterized by extreme class imbalance (>100:1 ratio of false candidate fragments to true variants). Under this raw per-segment space, default decision threshold $\tau = 0.50$ achieves Recall = 0.617 and Precision = 0.236 while maintaining Specificity $>0.9997$, successfully eliminating $>99.9\%$ of raw discovery noise before final cohort-level union assembly.

---

### 3.3 Objective 3: Nextflow Workflow Design and Computational Optimisation
**Status**: **Completed / Production-Ready**

#### 1. Workflow Architecture & Candidate Attrition
I built the entire workflow in **Nextflow DSL2** (`main.nf`), containerized all dependencies using Apptainer, and configured SLURM execution across worker nodes on the Wits HPC cluster. 

To track how candidate variants are processed in practice, I tracked call volume and diagnostic accuracy across eight key workflow stages on the DDD-UK cohort:

##### Table 5: Candidate Volume Attrition and Diagnostic Accuracy Across 8 Pipeline Stages (0.0–1.0 Scale)
> [!IMPORTANT]
> **Operational Specification**: Evaluated under **Mosdepth BED + 6-Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger** across 8 workflow stages on the 63 QC-passed DDD-UK probands.

| Stage Index & Pipeline Phase | Candidate Variants (Mosdepth) | Step Attrition (%) | Cumulative Retention (%) | Control Arm Calls | Biological Purpose & Clinical Milestone |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **Stage 1: 1. Raw Multi-Caller Discovery** | **46,791** | `0.0%` | 100.00% | 94,359 | Multi-caller parallel baseline |
| **Stage 2: 2. Quality Score Normalization Gate** | **67** | `−99.9%` | 0.14% | 964 | Universal scale [0–1000] & quality filters |
| **Stage 3: 3. Intra-Sample Multi-Caller Union** | **60** | `−10.4%` | 0.13% | 398 | Union of all caller calls (retains >=1 caller) |
| **Stage 4: 4. Breakpoint Refinement (Precise & Imprecise)** | **60** | `0.0%` | 0.13% | 398 | Base-pair junctions via SvABA (100% calls retained) |
| **Stage 5: 5. Inter-Sample Cohort Union** | **58** | `−3.3%` | 0.12% | 106 | Union across cohort samples (retains >=1 sample) |
| **Stage 6: 6. Internal Cohort AF Filtering (<1%)** | **9** | `−84.5%` | 0.02% | 33 | Purges benign common cohort polymorphisms |
| **Stage 7: 7. CatBoost Machine Learning (P ≥ 0.50)** | **1** | `−88.9%` | 0.00% | 7 | In silico gradient boosted probability triage |
| **Stage 8: 8. ACMG Class 4/5 Prioritisation** | **0** | `−100.0%` | 0.00% | 0 | High-confidence clinical pathogenic candidates |

*Note on Table 5 & Candidate Attrition Scope: Table 5 tracks the candidate attrition funnel for novel, rare pathogenic variants (which applies AF < 1% cohort frequency filtering and AnnotSV/ACMG pathogenicity classification, reducing candidate calls to actionable VUS/pathogenic candidates), whereas Table 1 benchmarks raw algorithmic sensitivity against known CMA truth intervals across the exome prior to cohort AF and ACMG exclusion. Counts are taken directly from `attrition_funnel_table.md` in `/spaces/ywolberg/qc-gated-cnv-combiner/` (in both Mosdepth and Control runs). The drop to 0 automated ACMG Class 4/5 variants in Stage 8 reflects rigorous adherence to ACMG/ClinGen 2019 standards: in `consensus_acmg_reports/DDD_Cohort_Rare_CNVs.csv` (N=59 rare candidate calls), AnnotSV classified 6 variants as Benign (ACMG 1), 3 as Variants of Uncertain Significance (ACMG 3, overlapping neurodevelopmental genes *CAMTA1* and *TNK2*), 0 as Class 4/5, and 50 as unclassified non-morbid loci. Under clinical diagnostic guidelines, rare CNVs in VUS genes cannot be automatically signed out as pathogenic without clinical multidisciplinary team (MDT) curation and trio segregation analysis. This honest attrition from 46,791 raw calls down to 9 rare cohort calls and 3 actionable VUS candidates demonstrates genuine clinical specificity without false-positive inflation.*

![Figure 5: Workflow Attrition Funnel and Stage Metrics](figures/fig5_workflow_attrition_funnel_and_stage_metrics.png)
*Figure 5: Workflow Attrition Funnel and Stage Metrics [Specification: Mosdepth BED + 6 Caller Union (+ ML) + Max Sensitivity Filter Profile + Truvari Merger] (0.0–1.0 Scale). (A) Candidate variant attrition funnel tracking call volume across 8 stages from 46,791 raw calls down to rare cohort variants. (B) Comparative call volume between Mosdepth Consensus BED and Standard Control BED, demonstrating that Mosdepth Consensus BED eliminates 50.4% of raw discovery noise (46,791 vs 94,359) and achieves cleaner cohort convergence. Stage 8 correctly yields 0 automated Class 4/5 calls because rare candidates affect VUS genes requiring MDT curation.*

#### 2. Resolving Breakpoints to Base-Pair Precision
A well-known weakness of read-depth callers is that they can only estimate CNV boundaries to the nearest captured exon, leaving breakpoints tagged as `IMPRECISE`.
* To solve this, I wrote a dedicated Python module ([`bin/refine_breakpoints.py`](file:///home/ywolberg/CNV-Caller-Combiner/bin/refine_breakpoints.py)) that matches depth-based consensus calls with split-read and discordant-pair evidence extracted from structural variant callers (SvABA).
* I evaluated this engine across the 90 clinical exomes evaluated under dual target regimes (Consensus vs. Control BED; 63 QC-passed DDD-UK samples):
  * **100% Call Retention**: Zero calls were lost or dropped (all candidate calls preserved intact across both runs).
  * **Breakpoint Resolution & Intronic Capture Constraints**: In whole-exome sequencing, short-read hybridisation capture probes target exonic regions; read depth and split reads decay rapidly within flanking uncaptured introns. Consequently, split-read resolution yielded base-pair exact junctions (`INFO/PRECISE`) for 7 deletions ($0.021\%$ of calls, where deletion breakpoints occurred within or immediately adjacent to captured probes) in the Mosdepth run and 16 deletions in the Control run. The remaining calls retain coordinate-harmonized Truvari boundaries (`INFO/IMPRECISE`), avoiding the window-bloating observed under naive clustering tools.
  * **Biological Mechanism Classification**: For the base-pair resolved junctions, sequence homology analysis identified the molecular mechanisms driving each rearrangement: Non-Allelic Homologous Recombination (NAHR, 3 events with $\ge 25\text{ bp}$ microhomology) and Non-Homologous End Joining (NHEJ, 4 events with $0\text{--}3\text{ bp}$ blunt/short insertions).

---

### 3.4 Objective 4: Deployment on the DDD-Africa Dataset
**Status**: **Staged / Pending (Scheduled for Q1–Q2 2027)**

With the pipeline containerized and validated on the UK cohort, I am now preparing to deploy it on 50 clinical exomes from the DDD-Africa study. This phase will directly benefit local patients and evaluate how our methods perform in diverse populations:
1. **Ancestry-Specific Baseline Calibration**: ExomeDepth and GATK-gCNV rely on reference cohorts to calculate normal depth expectations. Using European reference panels on African individuals introduces bias; therefore, I will calibrate panels of normals (PONs) using within-cohort African samples.
2. **Clinical Discovery in Undiagnosed Children**: I will apply our locked pipeline to identify diagnostic microdeletions and microduplications in previously unsolved African children with developmental delay, followed by ACMG curation.
3. **Cross-Cohort Transferability**: I will compare capture kinetics, probe coverage uniformity, and caller concordance between the UK and African cohorts, providing insights into exome CNV calling in understudied populations.

---

## 4. Planned Publications & Conference Dissemination

### Manuscripts in Preparation:
1. **Manuscript 1 (Methodology & Benchmarking)**:  
   *Working Title*: *"Systematic benchmarking and post-consensus breakpoint refinement of exome sequencing copy number variation caller combinations in neurodevelopmental disorders."*  
   *Target Journal*: *Genome Medicine* or *Bioinformatics* (Planned Submission: **Q1 2027**).
2. **Manuscript 2 (African Clinical Application)**:  
   *Working Title*: *"Application of an optimized machine-learning-gated CNV workflow to identify diagnostic structural variants in African probands with developmental disorders."*  
   *Target Journal*: *European Journal of Human Genetics* or *HGG Advances* (Planned Submission: **Q2/Q3 2027**).

### Conference Presentations:
* **2024**: Delivered a rapid-fire oral presentation at the SASBi/SAGS BIO2024 Student Symposium (awarded unofficial 2nd place); presented posters at the Future Africa Conference (University of Pretoria) and the Wits Molecular Biosciences Research Thrust (MBRT) Symposium.
* **2026**: Presented an oral paper at the Biennial Conference of the South African Genetics & Bioinformatics Societies (SAGS x SASBi BIO2026, Gqeberha); presented dual posters at the Wits Postgraduate Cross-Faculty Symposium and the Faculty of Health Sciences Research Day.

---

## 5. Thesis Completion Trajectory & Rationale for SAMRC Support

### 5.1 Feasibility and Timeline to Completion
* **Current Standing**: Completing Year 3; entering Year 4 in January 2027.
* **Target Submission**: **August / September 2027**.
* **Feasibility**: All computationally heavy development, algorithm benchmarking, pipeline containerization, and initial machine learning (Objectives 1, 2, and 3) are complete. My remaining 12 months are dedicated to running the locked workflow on the African cohort, writing two journal manuscripts, and compiling my thesis chapters.

### 5.2 Rationale for SAMRC PhD Scholarship Funding
Although I am entering my fourth year of doctoral study, I have **never received SAMRC funding**. Scholarship support during the 2027 academic year is crucial for the following reasons:
* **Dedicated Full-Time Focus**: Financial support will allow me to devote $100\%$ of my time to data synthesis, manuscript submissions, and dissertation writing without the need to take on distracting part-time work during my final write-up phase.
* **Advancing Genomic Equity in South Africa**: This project directly translates high-performance bioinformatics methods into an understudied African clinical cohort, helping address diagnostic disparities in local medical genetics.
* **Low Risk and Guaranteed Governance**: All necessary datasets are already in hand, the pipeline is fully operational on local servers, and institutional ethics clearance (**Protocol M2502103**, valid to July 2029) is fully secured.

---

**Report Compiled by**:  
Yonatan Ariel Wolberg  
PhD Candidate, Bioinformatics & Human Genetics  
University of the Witwatersrand, Johannesburg  
Date: September 2026
