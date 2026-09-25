# Breakpoint Refinement Module

## 1. Methodology & Algorithm
A fundamental limitation of read-depth based exome CNV callers is **breakpoint imprecision**: because callers infer copy number from discrete captured exons, the reported start and end coordinates correspond to the target exon boundaries rather than the true genomic breakends, frequently introducing thousands of base pairs of coordinate slack.

The **Breakpoint Refinement Module** (`bin/refine_breakpoints.py`, orchestrated by `REFINE_CONSENSUS_BREAKPOINTS`) bridges this gap by marrying depth-based consensus calls with nucleotide-level structural variant breakpoints derived from **SvABA** (or DRAGEN SV):
* **Size-Aware Search Windows**: Establishes adaptive flanking search margins around depth-based consensus calls. For focal variants (< 10 kb), smaller search margins are applied, while larger macro-events allow broader flanking windows.
* **Reciprocal Overlap Verification**: Requires a minimum reciprocal overlap (default: 20%) between the depth consensus event and the candidate split-read/assembly event to ensure they represent the same underlying biological lesion.
* **Coordinate Updating & Provenance Tracking**: When high-confidence split-read breakpoint evidence is identified, the VCF record's `POS` and `END` coordinates are updated to exact nucleotide positions. The INFO field is updated with `REFINED=TRUE`, `ORIGINAL_POS`, `ORIGINAL_END`, and `DELTA_BP`.

---

## 2. Nextflow Implementation & Processes
Encapsulated in `modules/refinement/breakpoint-refinement.nf`:

| Process | Container / Label | Function |
| :--- | :--- | :--- |
| `REFINE_CONSENSUS_BREAKPOINTS` | `pysam` | Matches depth consensus VCFs with sample-matched SvABA SV VCFs and updates coordinates using `bin/refine_breakpoints.py`. |

---

## 3. Inputs and References
* **Consensus VCF**: Intra-sample depth consensus VCF (`.vcf.gz` + `.tbi`) from SURVIVOR or Truvari.
* **SV VCF**: Sample-matched nucleotide-resolution SV VCF (`.vcf.gz` + `.tbi`) from SvABA (`*.svaba.sv.vcf.gz`) or DRAGEN SV.

---

## 4. Key Parameters & Thresholds
* `params.refinement_max_flank`: Maximum genomic distance for flanking search window (default: `50,000 bp`).
* `params.refinement_max_delta`: Maximum permitted shift between depth and assembly breakends (default: `25,000 bp`).
* `params.refinement_focal_threshold`: Size cutoff for focal CNVs (default: `10,000 bp`).
* `params.refinement_min_recip_overlap`: Minimum reciprocal overlap between depth call and SV event (default: `0.20`).

---

## 5. Output Files
* `consensus_refined/<SAMPLE_ID>.refined.vcf.gz`: Refined consensus VCF with updated nucleotide-resolution coordinates.
* `consensus_refined/<SAMPLE_ID>_refinement_summary.txt`: Diagnostic log detailing refined, unrefined, and ambiguous calls.
