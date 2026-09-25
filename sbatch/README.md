# Slurm Batch Submission Scripts (`sbatch/`)

This directory contains Slurm (`sbatch`) batch scripts tailored for submitting workflows and long-running pipeline tasks to the high-performance computing (HPC) cluster (specifically ZA-WITS-CORE).

---

## 📄 Core Submission Scripts

### 1. End-to-End Individual Callers
These scripts run individual CNV caller workflows (e.g., calling, filtering, ACMG, and evaluation) in end-to-end mode:
- **[run-workflow-canoes-e2e-90.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run-workflow-canoes-e2e-90.sbatch)**
- **[run-workflow-cnvkit-e2e-90.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run-workflow-cnvkit-e2e-90.sbatch)**
- **[run-workflow-gcnv-e2e-90.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run-workflow-gcnv-e2e-90.sbatch)**
- **[run-workflow-xhmm-e2e-90.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run-workflow-xhmm-e2e-90.sbatch)**
- *Variants ending in `-hard.sbatch` use stricter, hard-filtered sets.*

### 2. QC & General Pipelines
- **[run-qc-90.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run-qc-90.sbatch)**: Submits quality control/quality assessment steps.
- **[run_slurm_consensus_and_eval.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run_slurm_consensus_and_eval.sbatch)**: Runs the master consensus, merging, ML combiner, and final assessment workflow.

### 3. Utility & Maintenance
- **[move_results.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/move_results.sbatch)**: Maintenance script to clean up, sort, or copy pipeline run directories into downstream directories.
- **[run_combined_tests.sbatch](file:///home/ywolberg/CNV-Caller-Combiner/sbatch/run_combined_tests.sbatch)**: Automates running system integration and verification tests under Slurm allocation.

---

## ⚙️ How to Submit

To execute any script on the cluster, run:

```bash
sbatch sbatch/run_slurm_consensus_and_eval.sbatch
```
