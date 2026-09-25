# Execution & Testing Scripts (`scripts/`)

This directory contains utility shell scripts used to trigger local pipeline runs, batch testing, and failure debugging.

## Scripts

- `debug_last_failed.sh`: Inspects the most recent Nextflow `.exitcode` and `.command.err` in `.nextflow.log` / `work/` for fast troubleshooting.
- `run_local_combinations.sh`: Executes end-to-end combinations across multiple caller parameter configurations locally.
- `run_quick_test.sh`: Runs a lightweight fast validation on test samples to verify end-to-end functionality.

*(Backward compatibility symlinks are maintained at repository root).*
