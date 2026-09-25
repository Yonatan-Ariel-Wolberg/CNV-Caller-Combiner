# Truvari (Structural Variant Comparison and Benchmarking Engine)

## 1. Methodology & Algorithm
**Truvari** is a benchmark-grade structural variant comparison framework that incorporates geometric reciprocal overlap and sequence similarity scoring:
* **Reciprocal Overlap & Size Matching**: Evaluates genomic coordinate intersection using strict reciprocal overlap criteria (e.g., minimum 50% reciprocal overlap), ensuring that calls of vastly divergent sizes are not mistakenly lumped together.
* **`truvari collapse`**: Automatically collapses redundant or near-identical calls emitted by different algorithms into unified representative consensus variants.
* **Multi-Allelic Normalization**: Standardizes multi-allelic records using `bcftools norm -m -any` prior to collapse, preventing tool crashes and ensuring diploid representation.
* **Provenance & Metadata Aggregation**: Harmonizes INFO annotations across merged calls, propagating caller attribution (`CALLERS=CANOES,CLAMMS,XHMM`), caller count (`NUM_CALLERS`), and Phred-scaled consensus quality.

---

## 2. Nextflow Implementation & Processes
The Truvari module is encapsulated in `modules/sv-mergers/modules-truvari.nf` and includes the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `MERGE_VCFS` | `bcftools` | Concatenates sorted per-caller VCFs, normalizes sample identifiers in VCF headers, and splits multi-allelic records. |
| `TRUVARI_COLLAPSE` | `truvari` | Executes `truvari collapse` using specified reciprocal overlap and distance thresholds. |
| `ANNOTATE_TRUVARI_CONSENSUS` | `python` | Annotates output VCFs with caller provenance tags, supporting counts, and consensus confidence metrics. |
| `BGZIP_SORT_INDEX_VCF` | `bcftools` | Compresses and creates tabix indices for Truvari consensus outputs. |

---

## 3. Inputs and Reference Requirements
* **Input Callsets**: Filtered, coordinate-sorted VCFs from individual callers per sample.
* **Reference Genome**: GRCh38/hg38 FASTA file and `.fai` index.

---

## 4. Key Parameters & Thresholds
* `params.truvari_pctovl`: Minimum reciprocal overlap percentage required to collapse calls (default: `0.50`).
* `params.truvari_refdist`: Maximum distance allowed between reference coordinates (default: `500 bp`).
* `params.truvari_sizemin`: Minimum variant size threshold (default: `50 bp`).

---

## 5. Output Files
* `consensus_truvari/<SAMPLE_ID>_truvari_collapsed.vcf.gz`: Consolidated consensus VCF.
* `consensus_truvari/<SAMPLE_ID>_truvari_bench.json`: Truvari collapse summary statistics and overlap logs.
