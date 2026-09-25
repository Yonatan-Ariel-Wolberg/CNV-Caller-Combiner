# Workflow Utilities (`workflows/utilities/`)

This directory contains specialized auxiliary Nextflow workflows and utility pipelines used for cohort preprocessing, mosdepth coverage calculations, score-separated evaluations, and cloud platform setups.

## Included Workflows

- `ddduk-90-control-filter-merge-freq.nf`: Filters and merges allele frequencies across control cohort panels.
- `ddduk-90-mosdepth-filter-merge-freq.nf`: Computes target coverage with mosdepth, applies threshold filters, and computes frequencies.
- `evaluate_score_separated.nf`: Evaluates caller and combiner performance stratified by score bands and thresholds.
- `icav2_config_set.nf`: Generates and verifies Illumina Connected Analytics (ICAv2) configuration and project tokens.

*(Backward compatibility symlinks are maintained at repository root).*
