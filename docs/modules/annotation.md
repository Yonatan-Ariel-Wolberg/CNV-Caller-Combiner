# Clinical Annotation & ACMG Stratification Module (AnnotSV)

## 1. Methodology & Algorithm
Clinical translation of CNV findings requires rigorous interpretation adhering to established genetic standards. The **Annotation Module** integrates **AnnotSV**, the premier structural variant annotation engine, to classify variants according to the **ACMG / ClinGen 2019/2020 Standards and Guidelines for Germline CNVs**:
* **Joint Gene-Level & Variant-Level Annotation**: Evaluates each CNV across both its full genomic span ("split" mode) and on an overlapping gene-by-gene basis ("full" mode), annotating impacted transcripts, coding exons, promoters, and untranslated regions (UTRs).
* **Decoupled Architecture**: Automatically downloads, stages, and caches the multi-gigabyte human annotation databases locally, verifying archive checksums and preventing redundant network downloads across Nextflow resumes.
* **ACMG Pathogenicity Scoring**: Computes evidence-based scores for deletions (haploinsufficiency, HI) and duplications (triplosensitivity, TS) across 5 standard clinical tiers:
  - **Class 5**: Pathogenic (high penetrance, established dosage sensitivity, microdeletion/microduplication syndromes).
  - **Class 4**: Likely Pathogenic.
  - **Class 3**: Variant of Uncertain Significance (VUS).
  - **Class 2**: Likely Benign.
  - **Class 1**: Benign.
* **Clinical Prioritization Stratification**: Generates stratified VCF files and spreadsheet-ready TSV tables (`*ACMG_4_5.tsv`) that isolate actionable Pathogenic and Likely Pathogenic variants for medical geneticists.

---

## 2. Nextflow Implementation & Processes
Encapsulated in `modules/annotation/annotation.nf`:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `DOWNLOAD_ANNOTSV_DB` | `annotsv` | Clones AnnotSV repository, downloads annotation databases in parallel, verifies archive integrity, and compiles local DB cache. |
| `RUN_ANNOTSV` | `annotsv` | Annotates filtered cohort VCFs against the local database, outputting full tabular annotation tables (`*.annotated.tsv`). |
| `ADD_ACMG_ANNOTATIONS` | `python` | Translates AnnotSV annotations into VCF INFO fields (`ACMG_CLASS`, `GENES`, `HI_SCORE`, `TS_SCORE`). |
| `STRATIFY_VCF_BY_ACMG` | `bcftools` | Splits the annotated cohort VCF into discrete clinical tier subsets (`ACMG_4_5`, `ACMG_5`, `ACMG_4`, `ACMG_3`, `ACMG_1_2`). |

---

## 3. Output Files
* `consensus_cohort/cohort_rare.annotated.tsv`: Complete multi-column AnnotSV annotation table.
* `consensus_evaluate/<RUN_NAME>_ACMG_4_5.tsv`: Spreadsheet-ready high-confidence clinical candidate list.
* `consensus_evaluate/<RUN_NAME>_ACMG_5.vcf.gz`: Dedicated Pathogenic tier VCF.
* `consensus_evaluate/<RUN_NAME>_ACMG_4.vcf.gz`: Dedicated Likely Pathogenic tier VCF.
