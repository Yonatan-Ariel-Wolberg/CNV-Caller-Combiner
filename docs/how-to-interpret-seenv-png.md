# How to Interpret SeeNV CNV Infographics (`.png` Cards)

The **CNV-Caller-Combiner** pipeline generates publication-ready, standardized visual infographic cards for candidate CNVs using either the official **SeeNV** binary (`seenv -p`, Bradshaw et al., *J Mol Diagn* 2023) or the native high-performance cohort-contrast visualizer ([bin/render_cnv_infographics.py](file:///home/ywolberg/CNV-Caller-Combiner/bin/render_cnv_infographics.py)).

These cards provide rapid, at-a-glance orthogonal validation of candidate copy number variants by contrasting the target sample's exonic read depth against the **Cohort Reference Panel (Option 2)**.

**Why the Cohort Reference Panel (Option 2)?**:
Rather than relying on generic external reference panels (e.g. 1000 Genomes WGS), our pipeline constructs a baseline directly from the cohort's own sequenced BAMs. This perfectly normalizes for capture kit hybridisation kinetics, GC-content biases, sequencer batch effects, and target library prep conditions.

---

## 1. Infographic Card Architecture

Each generated `.png` card is partitioned into three coordinated evidence panels:

```text
+-----------------------------------------------------------------------------------------+
|                                    SeeNV CNV CARD                                       |
|  SAMPLE: <SAMPLE_ID>      |  VARIANT: DEL (45,210 bp)  |  LOCUS: chr2:15,200,100-15,245,310|
+---------------------------------------------------+-------------------------------------+
|                                                   | PANEL 2: CLINICAL EVIDENCE CARD     |
| PANEL 1: RELATIVE EXONIC DEPTH PROFILE            | Clinical Curation Evidence:         |
|                                                   |   Sample ID:        <SAMPLE_ID>         |
| Ratio                                             |   Variant:          DEL (45,210 bp)     |
| 2.0 |                                             |   Coordinates:      chr2:15.2M-15.24M   |
|     |                                             |   Consensus:        CANOES,XHMM,GATK    |
| 1.5 |                                             |   Supporting Tools: 3 callers           |
|     |                                             |   Depth Ratio:      0.48x (Observed)    |
| 1.0 | ===== [ 1.02x ] ============= [ 0.98x ] === |                                     |
|     |  (Upstream Flank)             (Down Flank)  |   [ACMG Class 5: Pathogenic]            |
| 0.5 | .................. [ 0.48x ] .............. |-------------------------------------|
|     |                   (Target DEL)              | PANEL 3: TECHNICAL ASSESSMENT       |
| 0.0 +-------------------------------------------+ | ✓ TECHNICAL EVIDENCE: SUPPORTED     |
|        Flank Left     CNV Target     Flank Right  | • Depth drop concordant with DEL    |
|                                                   | • Baseline Flank Depth: 52.4x       |
|                                                   | • Target Region Depth:  25.1x       |
|                                                   | • Action: Proceed to sign-off       |
+---------------------------------------------------+-------------------------------------+
```

---

## 2. Detailed Panel Interpretation

### Panel 1: Relative Exonic Depth Profile (Left Main Chart)
This panel displays the quantitative read-depth signature across three contiguous windows:
1. **Upstream Flank (Left Bar)**: The un-mutated diploid genomic region immediately preceding the variant.
2. **CNV Target (Center Bar)**: The precise interval of the predicted deletion or duplication.
3. **Downstream Flank (Right Bar)**: The un-mutated diploid genomic region immediately following the variant.

#### Reference Baseline Lines
- **Diploid Baseline (1.0x, Grey Dashed Line)**: Represents expected normal diploid copy number ($2n$). Flanking regions should hover closely around $1.0\times$ (typically $0.85\times - 1.15\times$).
- **Expected Ratio Line (Dotted Line)**:
  - **Heterozygous Deletion ($1n$)**: Red dotted line at **$0.50\times$**.
  - **Heterozygous Duplication ($3n$)**: Blue dotted line at **$1.50\times$**.
  - **Homozygous Deletion ($0n$)**: Baseline at **$0.00\times$**.

#### Bar Labels
- Top number indicates the **Normalized Depth Ratio** ($\text{Target Coverage} / \text{Flanking Baseline}$).
- Bottom number indicates the **Raw Mean Coverage** (e.g. `25.1x cov`).

---

### Panel 2: Clinical Curation & Evidence Metadata (Top-Right Card)
This metadata card provides instant clinical context without requiring navigation across multiple VCF files:
- **Sample ID**: Target sample accession.
- **Variant & Size**: CNV type (DEL/DUP) and length in base pairs.
- **Consensus Callers**: Names of all orthogonal algorithms that called this event (e.g. `CANOES,CLAMMS,XHMM,GATK_GCNV,CNVKIT,EXOMEDEPTH`).
- **Supporting Tools**: Total count of concurring callers. Higher counts significantly correlate with true positive status.
- **Depth Ratio**: Raw observed target-to-baseline ratio.
- **ACMG Tier Badge**:
  - `[ACMG Class 5: Pathogenic]` (Dark Red): High clinical priority; overlaps established haploinsufficient or triplosensitive disease gene.
  - `[ACMG Class 4: Likely Pathogenic]` (Orange): Clinically significant event with substantial disease overlap.
  - `[ACMG Class 3: VUS]` / `[ACMG Evaluated]` (Blue): Variant of uncertain significance.

---

### Panel 3: Technical Validation Recommendation (Bottom-Right Card)
Automated heuristic engine evaluating whether the physical read depth matches the algorithmic call:
- **Status Banner**:
  - `✓ TECHNICAL EVIDENCE: SUPPORTED` (Green): Target depth ratio satisfies calibrated thresholds ($\le 0.65\times$ for DEL, $\ge 1.30\times$ for DUP).
  - `⚠ TECHNICAL EVIDENCE: CHECK FLANKS` (Orange/Red): Ratio is ambiguous, flanking coverage is low ($< 15\times$), or depth ratio conflicts with the caller label.
- **Flanking vs Target Depth**: Displays raw un-normalized coverage across both regions to confirm adequate sequencing depth.
- **Recommended Action**: Direct instructions for the analyst:
  - *"Proceed to clinical sign-off"* $\rightarrow$ Call is robust, clear, and unambiguous.
  - *"Inspect in IGV for split reads"* $\rightarrow$ Call has borderline coverage or single-tool support; requires base-level verification in the IGV HTML report.

---

## 3. Cohort Dispersion Plots (Official SeeNV Mode)

When official SeeNV is executed with the pre-computed cohort reference database (`SEENV_BUILD_COHORT_DB`), the left panel renders an **Exon-by-Exon Multi-Sample Scatter Plot**:

```text
Exon Depth Ratio
1.5 |    o       o            o       o   (Duplication Level)
    |  +---+   +---+        +---+   +---+
1.0 | -|   |---|   |--------|   |---|   |- (Diploid Cohort Envelope: 5th-95th percentile)
    |  +---+   +---+  +---+ +---+   +---+
0.5 |                 |   |               <-- Target Sample Dots drop to ~0.5x
    |                 +---+
0.0 +------------------------------------
       Exon 1  Exon 2 Exon 3 Exon 4 Exon 5
```

- **Grey Box / Shaded Envelope**: The distribution of normalized read depth across all **90 cohort samples** at each individual exon.
- **Solid Colored Points**: The target sample's coverage at that specific exon.
- **True Positive Criteria**:
  - All exons within the CNV interval drop/rise simultaneously outside the cohort normal envelope.
  - Flanking exons of the same sample remain centered inside the cohort interquartile range.

---

## 4. Clinical Decision & Interpretation Matrix

| Metric | True Heterozygous Deletion | True Heterozygous Duplication | Common Artifact / Noise |
| :--- | :--- | :--- | :--- |
| **Normalized Ratio** | $0.35\times - 0.65\times$ | $1.30\times - 1.70\times$ | $0.75\times - 1.20\times$ |
| **Flanking Depth** | Adequate ($\ge 20\times$) and stable | Adequate ($\ge 20\times$) and stable | Very low ($< 10\times$) or noisy |
| **Cohort Dispersion** | Sample is a clear outlier below cohort 5th percentile | Sample is a clear outlier above cohort 95th percentile | Sample falls inside normal cohort variation band |
| **Caller Support** | $\ge 2$ callers (e.g. CANOES + XHMM) | $\ge 2$ callers (e.g. CNVkit + GATK) | 1 single caller |
| **ACMG Priority** | Class 4 or 5 | Class 4 or 5 | Class 1, 2, or non-coding |

---

## 5. Summary: Synergistic Use of SeeNV PNG and IGV HTML

| Task | Primary Tool | Why? |
| :--- | :--- | :--- |
| **Rapid Cohort Screening** | **SeeNV PNG** | Instant visual comparison against all 90 cohort samples to rule out cohort-wide batch effects in seconds. |
| **Clinical Report Figures** | **SeeNV PNG** | Clean, self-contained, publication-ready vector image for multidisciplinary team (MDT) meetings. |
| **Single-Nucleotide Validation** | **IGV HTML** | Inspect split reads (soft-clips) and discordant read pairs at the exact breakpoint junctions. |
| **Pseudogene / Repeat Checking** | **IGV HTML** | Confirm mapping quality ($\text{MAPQ} = 60$ vs $0$) and rule out homologous gene cross-mapping. |
