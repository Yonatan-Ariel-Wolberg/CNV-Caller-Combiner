# Illumina DRAGEN (Dynamic Read Analysis for GENomics) via ICAv2

## 1. Methodology & Algorithm
**Illumina DRAGEN** is a high-throughput, hardware-accelerated bioinformatics pipeline platform offering high-accuracy variant discovery. In the **CNV-Caller-Combiner** pipeline, DRAGEN is integrated through an automated cloud execution wrapper using the **Illumina Connected Analytics version 2 (ICAv2)** CLI and API:
* **Target-Aware Read-Depth Model**: Normalizes coverage per targeted region against GC content, regional mappability, and cohort reference baselines.
* **Integrated Structural Variant (SV) Detection**: Concurrently leverages paired-end insert size distributions and split-read mapping evidence to detect deletions, duplications, and complex rearrangements.
* **Cloud API Orchestration**: Automatically stages local BAM/CRAM files into secure cloud project storage, launches DRAGEN Germline Enrichment pipelines on FPGA hardware, monitors job completion status, and streams back final VCF outputs.

---

## 2. Nextflow Implementation & Processes
The DRAGEN ICAv2 module is encapsulated in `modules/callers/modules-icav2-dragen.nf` and includes the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `ICAV2_UPLOAD_INPUTS` | `icav2_cli` | Verifies and uploads alignment files (`.bam`/`.cram`) and index files (`.bai`/`.crai`) into the designated ICAv2 cloud storage project. |
| `ICAV2_LAUNCH_DRAGEN` | `icav2_cli` | Submits the DRAGEN Germline Enrichment analysis workflow with targeted intervals and reference specifications. |
| `ICAV2_POLL_STATUS` | `icav2_cli` | Monitors the asynchronous cloud execution until completion (`SUCCEEDED`, `FAILED`, or `ABORTED`). |
| `ICAV2_DOWNLOAD_RESULTS` | `icav2_cli` | Streams back output VCF files (`.cnv.vcf.gz`, `.sv.vcf.gz`), metrics, and coverage summaries from cloud storage to local results folders. |
| `FORMAT_DRAGEN_VCF` | `python` | Standardizes DRAGEN CNV format attributes, ensuring compatibility with downstream SURVIVOR and Truvari mergers. |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM or CRAM files.
* **Target Intervals**: Target BED file defining captured exons.
* **ICAv2 Credentials**: `ICAV2_API_KEY` exported in the execution environment.
* **Cloud Project**: Pre-configured ICAv2 project ID and storage container.

---

## 4. Key Parameters & Quality Thresholds
* `params.icav2_project_id`: Cloud project identifier.
* `params.filter_dragen_qscore`: Quality score threshold (default: `>= 30`).
* `params.filter_dragen_min_size`: Minimum CNV event length in base pairs (default: `1000 bp`).

---

## 5. Output Files
* `out_DRAGEN/raw/<SAMPLE_ID>.dragen.raw.vcf.gz`: Raw output CNV VCF downloaded from ICAv2.
* `out_DRAGEN/filtered/<SAMPLE_ID>.dragen.filtered.vcf.gz`: Quality-filtered and normalized callset.
