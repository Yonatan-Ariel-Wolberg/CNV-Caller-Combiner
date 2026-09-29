# How to Interpret Headless IGV HTML Reports (`igv-reports`)

The **CNV-Caller-Combiner** pipeline generates standalone, zero-retention interactive HTML reports using `igv-reports` (container: `quay.io/biocontainers/igv-reports:1.16.0--pyh7e72e81_0`).

These reports embed micro-slices of alignment data (BAM reads), exome capture target intervals (probes BED), gene annotations (RefSeq/GENCODE), and clinical variant predictions into a single self-contained `.html` file.

**Zero Data Leakage Guarantee**: The IGV HTML report runs **100% offline in your local web browser**. It never uploads, transfers, or pings any genomic data or metadata to external servers or cloud services. You can safely open and inspect sensitive patient datasets on an air-gapped clinical workstation.

---

## 1. User Interface Overview

When you open the generated `.html` file in any modern web browser (Chrome, Firefox, Safari, Edge), you will see two primary components:

```text
+-----------------------------------------------------------------------------------------+
| [Search...] [Filter by Class: ACMG 5/4/3] [Export CSV]                                   |
| Table: SampleID | Chrom | Start | End | Type | Size | Callers | ACMG Class | Gene       |
| > <SAMPLE_ID> | chr17 | 43044295 | 43170245 | DEL | 125,950 bp | 5 | Class 5 | BRCA1        |
+-----------------------------------------------------------------------------------------+
|                                    IGV EMBEDDED VIEWER                                  |
| [Locus View: chr17:43,000,000-43,220,000]           [Breakpoint View: chr17:43,044,295]  |
| --------------------------------------------------------------------------------------- |
| Target Capture Probes Track:  |==|    |==|    |==|            |==|    |==|              |
| Coverage Histogram Track:     ████    ████    ████            ████    ████              |
| BAM Alignment Track:          ===>    ===>    ===>            ===>    ===>              |
| Gene Models Track:            BRCA1 ------> [Exon 1] -------- [Exon 2] -------->        |
+-----------------------------------------------------------------------------------------+
```

1. **Variant Navigation Table (Top Panel)**:
   - Lists every prioritized candidate CNV evaluated in the run.
   - Columns include: `CHROM`, `START`, `END`, `SVTYPE` (DEL or DUP), `SIZE_BP`, `CALLERS` (e.g. `CANOES,CLAMMS,XHMM,GATK`), `NUM_CALLERS`, `ACMG_CLASS`, and affected gene symbols.
   - **Clicking any row** automatically navigates the embedded IGV viewer below directly to that locus.
   - Search bar allows instant filtering by gene name (e.g., `DMD`, `BRCA1`), chromosome, or variant type.

2. **Dual Genomic Viewers (Bottom Panel)**:
   - **Wide Locus View (Left)**: Displays the entire CNV plus 25–50 kb of flanking genomic sequence to evaluate baseline coverage contrast.
   - **Breakpoint Detail View (Right)**: Zooms into the exact predicted breakpoint boundaries to search for split reads, microhomology, and discordant mate pairs.

---

## 2. Track-by-Track Interpretation

### Track A: Exome Target Probes (`probes.bed`)
- **Visual Representation**: Horizontal colored bars marking regions targeted by the exome capture kit (e.g., Twist, Agilent SureSelect, Illumina Exome).
- **Interpretation Rule**:
  - In Whole Exome Sequencing (WES), read coverage exists **only** where capture probes hybridize.
  - **Critical Validation Check**: If coverage drops in a region that has **no capture probes**, that drop is an expected **unprobed intronic/intergenic gap**, *not* a genomic deletion.
  - A genuine CNV must exhibit depth changes across **probed target intervals**.

---

### Track B: Depth of Coverage Histogram
- **Visual Representation**: Continuous grey/blue area histogram depicting local read depth at every nucleotide.
- **Copy Number Evaluation**:
  - **Diploid Baseline (Normal $2n$)**: Consistent with cohort average (typically $30\times - 100\times$).
  - **Heterozygous Deletion ($1n$)**: Mean depth drops to approximately **$0.40\times - 0.60\times$** of the flanking probed exons.
  - **Homozygous Deletion ($0n$)**: Read depth drops to **$0\times$** across all probed exons within the interval (with negligible background noise).
  - **Heterozygous Duplication ($3n$)**: Read depth increases to approximately **$1.35\times - 1.65\times$** relative to flanking probed exons.
  - **High-Copy Amplification ($\ge 4n$)**: Read depth exceeds **$\ge 2.0\times$** baseline.

```text
Expected Exonic Depth Profile:
Diploid (2n):      [ 60x ] ---- [ 62x ] ---- [ 58x ] ---- [ 61x ]
Het Deletion (1n): [ 60x ] ---- [ 28x ] ---- [ 31x ] ---- [ 62x ]  <-- ~50% drop
Het Duplication:   [ 60x ] ---- [ 92x ] ---- [ 89x ] ---- [ 61x ]  <-- ~50% increase
```

---

### Track C: BAM Alignment Track (Individual Reads)
Individual reads provide the nucleotide- and molecule-level proof supporting or refuting a read-depth call.

#### 1. Insert Size Coloring (Discordant Read Pairs)
In paired-end sequencing, fragments typically span 300–500 bp. IGV colors reads with abnormal inferred insert sizes:
- **Red Reads (Larger Insert Size than Expected)**:
  - **Interpretation: Structural Deletion**.
  - When a deletion occurs in the patient's DNA, the sequencer reads the flanking ends of the fragment. When aligned back to the reference genome, the two paired mates align abnormally far apart, triggering the red color.
- **Blue Reads (Smaller Insert Size than Expected)**:
  - **Interpretation: Structural Insertion or Short Fragment**.
- **Dark Blue / Cyan / Teal Reads (Inverted Orientation)**:
  - Expected orientation is Forward-Reverse (`FR`, $\rightarrow \leftarrow$).
  - Orientation `FF` ($\rightarrow \rightarrow$), `RR` ($\leftarrow \leftarrow$), or `RF` ($\leftarrow \rightarrow$) indicates an **Inversion** or a **Tandem Duplication breakpoint junction**.

#### 2. Soft-Clipping and Split Reads
- **Visual Representation**: Reads with multi-colored vertical stripes at their extremities.
- **Mechanism**: A split read spans the exact breakpoint junction. One portion of the read aligns to the start of the CNV, while the remaining portion belongs to the distal junction.
- **Clinical Significance**:
  - Soft-clipped reads aligning at the exact same base on both sides of a predicted CNV provide **single-nucleotide resolution proof** of the structural variant.

#### 3. Mapping Quality (MAPQ)
- **High MAPQ (MAPQ = 60)**: Opaque, solid reads. Aligns uniquely and unambiguously to the reference.
- **Low / Zero MAPQ (MAPQ = 0)**: Faded, transparent, or white reads. Multi-mapping reads that align equally well to multiple genomic locations.
- **Interpretation Rule**: If a suspected deletion or duplication is supported only by MAPQ 0 reads in a segmental duplication or pseudogene region (e.g., *SMN1/SMN2*, *GBA/GBAP1*, *CYP21A2*), exercise extreme caution: it may represent mapping artifacts rather than a true copy number change.

#### 4. B-Allele Frequency (BAF) of Heterozygous SNPs
- Inspect any single-nucleotide variants inside the CNV window:
  - **Heterozygous Deletion ($1n$)**: Complete **Loss of Heterozygosity (LOH)**. All informative SNPs will appear strictly homozygous ($0\%$ or $100\%$, reference or alternate). A true heterozygous deletion cannot contain genuine heterozygous SNPs ($50/50$ allele fraction).
  - **Heterozygous Duplication ($3n$)**: Allelic imbalance shifts heterozygous SNPs from $50:50$ to **$33:67$ or $67:33$** (e.g., $AAB$ or $ABB$).

---

### Track D: Gene Annotations & Transcripts
- Displays transcript models (exons as solid blocks, introns as lines with directional arrowheads indicating $5' \rightarrow 3'$ transcription).
- **Interpretation Rule**:
  - Determine whether the CNV covers the **whole gene**, **multiple contiguous exons**, or is **intragenic** (e.g. single-exon deletion).
  - Check whether the CNV disrupts the open reading frame (frameshift vs in-frame deletion) and whether it affects canonical transcripts.

---

## 3. Clinical Curation Decision Matrix

Use the following step-by-step checklist to review and sign off on variants:

| Verification Step | High-Confidence True Positive | Suspicious / Likely Artifact |
| :--- | :--- | :--- |
| **1. Probes Check** | All affected exons contain verified target capture probes. | Coverage drop occurs in an unprobed region or intronic gap. |
| **2. Coverage Step** | Sharp, uniform $\sim 50\%$ drop (DEL) or $\sim 50\%$ rise (DUP) bounded by probe coordinates. | Gradual sloping coverage or noisy fluctuations seen in all cohort samples. |
| **3. Breakpoint Junctions** | Split reads (soft-clips) and discordant read pairs (red/RF) align at the predicted ends. | No discordant or split reads present anywhere near the predicted boundaries. |
| **4. Mapping Quality** | Reads show high mapping quality ($\text{MAPQ} \ge 30$). | Reads are transparent/white ($\text{MAPQ} = 0$) overlapping known pseudogenes or repeats. |
| **5. Allele Balance** | Complete absence of heterozygous SNPs inside deletion (LOH observed). | Multiple high-confidence $50/50$ heterozygous SNPs present within the deleted interval. |
| **6. Caller Consensus** | $\ge 2$ orthogonal callers agree (e.g. CANOES + XHMM + GATK). | Single caller call without corroboration from any other tool. |

---

## 4. Common Edge Cases and Pitfalls

1. **Single-Exon Deletions with No Split Reads**:
   - In exome sequencing, breakpoints frequently fall in deep intronic regions not captured by the probes. A true single-exon deletion may therefore lack split reads in the BAM.
   - **Resolution**: Cross-reference with the [SeeNV Infographic](file:///home/ywolberg/CNV-Caller-Combiner/docs/how-to-interpret-seenv-png.md) to confirm whether the exon's normalized read depth is significantly depleted compared to the remaining 89 cohort samples.

2. **Extreme GC Content**:
   - Exons with $>75\%$ GC or $<25\%$ GC often suffer from PCR amplification dropouts.
   - Check if adjacent unrelated samples also show low coverage at that specific exon. If all samples drop, it is a capture/sequencing artifact.

3. **Duplication Orientation**:
   - Tandem duplications in direct head-to-tail orientation produce read pairs with inverted orientation (`RF`, mapped outwardly).
   - Interspersed duplications (translocated insertions) will show read pairs pointing to completely different chromosomes.
