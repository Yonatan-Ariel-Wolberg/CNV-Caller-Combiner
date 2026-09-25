# Parameter Configurations (`params/`)

This directory contains Nextflow parameter configuration files formatted in JSON, specifying paths, reference genomes, target BED intervals, output directories, and filtering thresholds across callers and workflows.

---

## 📁 Structure

### 1. Generic Parameter Templates (`params/general/`)
Reusable parameter templates that serve as blueprints for custom cohorts:
- **`params-full.json`**: End-to-end multi-caller and consensus execution parameters.
- **`params-survivor.json`** / **`params-truvari.json`**: Multi-caller intra-sample and inter-sample consensus merging parameters.
- **`params-filter-balanced.json`**: Balanced quality and size filtering profile (recommended default, optimizes F1-score).
- **`params-filter-high-precision.json`**: Stringent quality filtering profile (minimizes false positives, ideal for qPCR validation).
- **`params-filter-max-sensitivity.json`**: Relaxed quality filtering profile (maximizes variant discovery, ideal for undiagnosed rare disease cohorts).
- **Individual Caller Templates**:
  - `params-canoes.json`: CANOES negative binomial HMM calling parameters.
  - `params-clamms.json`: CLAMMS depth and mappability calling parameters.
  - `params-cnvkit.json`: CNVkit on/off-target segmentation parameters.
  - `params-exomedepth.json`: ExomeDepth beta-binomial reference selection parameters.
  - `params-gatk-gcnv.json`: GATK gCNV Bayesian interval model parameters.
  - `params-icav2-dragen.json`: Illumina DRAGEN cloud ICAv2 pipeline parameters.
  - `params-xhmm.json`: XHMM PCA normalization and HMM parameters.

### 2. Tool-Specific Runtime Configurations
- **[xhmm.conf](file:///home/ywolberg/CNV-Caller-Combiner/params/xhmm.conf)**: Runtime threshold definitions, PCA standard deviation variance cutoffs, and matrix normalization filters used by the XHMM engine.

---

## 🔧 Usage

When executing Nextflow pipelines, pass any parameter configuration file using `-params-file`:

```bash
# Run full consensus pipeline with balanced filtering
nextflow run cnv-caller-combiner.nf \
  -params-file params/general/params-full.json \
  -profile wits_cluster \
  --workflow full \
  --ml_pathway A \
  -resume
```

