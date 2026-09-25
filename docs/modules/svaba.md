# SvABA (Structural Variation and Indel Analysis by Local Assembly)

## 1. Methodology & Algorithm
**SvABA** is a specialized structural variant (SV) and indel caller that operates through local de novo sequence assembly rather than relying solely on read-depth fluctuations:
* **Local String-Graph Assembly**: Extracts discordant paired-end reads and split reads mapping around targeted regions. It constructs local string graphs to assemble contigs de novo without requiring prior knowledge of breakpoint positions.
* **Breakpoint Realignment with BWA-MEM**: Assembled contigs are realigned to the reference genome using an integrated BWA-MEM engine. Breakpoints are pinpointed at **single-nucleotide resolution**, revealing microhomology, non-templated sequence insertions, and exact genomic coordinates.
* **Target-Restricted Execution for Exomes**: In the **CNV-Caller-Combiner** pipeline, SvABA is constrained to the exome target BED intervals (`-k target_bed`) to optimize computational runtime while capturing high-confidence breakpoints across all targeted genes.
* **Breakpoint Refinement Engine**: SvABA SV calls serve as the primary nucleotide-resolution evidence source in the `REFINE_CONSENSUS_BREAKPOINTS` module, updating the approximate window boundaries of read-depth callers to precise genomic coordinates.

---

## 2. Nextflow Implementation & Processes
The SvABA module is encapsulated in `modules/callers/modules-svaba.nf` and includes the following atomic processes:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `RUN_SVABA_ASSEMBLY` | `svaba` | Prepares alignment inputs (handling BAM vs CRAM target slicing), ensures BWA index presence, and executes `svaba run` in target-restricted germline mode. |
| `BGZIP_SORT_INDEX_SV_TUPLE_VCF` | `bcftools` | Sorts, filters for canonical chromosomes, bgzips, and creates `.tbi` indices for raw SvABA structural variant VCFs (`*.svaba.sv.vcf.gz`). |

---

## 3. Inputs and Reference Requirements
* **Alignments**: Coordinate-sorted, indexed BAM or CRAM files (`.bam` + `.bai` or `.cram` + `.crai`).
* **Target Intervals**: Target BED file defining captured exons.
* **Reference Genome**: GRCh38/hg38 FASTA file with `.fai` index and BWA index files (`.bwt`, `.pac`, `.ann`, `.amb`, `.sa`). If BWA indices are missing, the process builds them automatically.

---

## 4. Key Parameters & Thresholds
* `params.svaba_target_padding`: Flanking margin added around target intervals (default: `100 bp`).
* `params.svaba_min_mapq`: Minimum mapping quality for assembly reads (default: `20`).

---

## 5. Output Files
* `out_SVABA/<SAMPLE_ID>.svaba.sv.vcf.gz`: Nucleotide-resolution structural variant calls with exact breakends (`BND`), deletions (`DEL`), and duplications (`DUP`).
* `out_SVABA/<SAMPLE_ID>.svaba.indel.vcf.gz`: Short indel calls derived from local assembly.
