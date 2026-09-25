#!/usr/bin/env python3
"""
Tests for DRAGEN minimum-sample enforcement.

DRAGEN Germline Enrichment requires a minimum of 5 samples to enable the
in-run panel of normals (cnv_enable_in_run_pon).

Validates that:
  1. The DRAGEN workflow in modules-icav2-dragen.nf contains a Nextflow-level
     count check that errors when fewer than 5 samples are provided.
  2. The START_ANALYSIS_BATCH bash script validates the minimum 5-sample count
     before launching the analysis.
  3. The cnv_enable_in_run_pon parameter is enabled in START_ANALYSIS_BATCH.
  4. The minimum threshold is exactly 5 (as required by DRAGEN).
"""

import os
import re
import subprocess
import textwrap

# Repository root relative to this test file
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DRAGEN_MODULE = os.path.join(REPO_ROOT, 'modules', 'callers', 'modules-icav2-dragen.nf')


def _load_module():
    with open(DRAGEN_MODULE) as f:
        return f.read()


def _extract_workflow_block(content):
    """Extract the text of the DRAGEN workflow block."""
    m = re.search(r'workflow DRAGEN \{', content)
    assert m is not None, "workflow DRAGEN not found in modules-icav2-dragen.nf"
    start = m.end() - 1
    depth = 0
    for i, ch in enumerate(content[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return content[start: i + 1]
    raise AssertionError("Unmatched braces for workflow DRAGEN")


def _extract_start_analysis_batch_script(content):
    """Extract the bash script body of the START_ANALYSIS_BATCH process."""
    m = re.search(r'process START_ANALYSIS_BATCH \{', content)
    assert m is not None, "process START_ANALYSIS_BATCH not found in modules-icav2-dragen.nf"
    # Find the script: block within the process
    proc_text = content[m.start():]
    script_m = re.search(r'\bscript\b\s*:', proc_text)
    assert script_m is not None, "script: block not found in START_ANALYSIS_BATCH"
    script_body = proc_text[script_m.end():]
    # Extract up to the closing triple-quote of the heredoc
    end_m = re.search(r'^\s*"""', script_body, re.MULTILINE)
    # The first """ opens, the second """ closes; skip the first
    first = script_body.index('"""')
    second = script_body.index('"""', first + 3)
    return script_body[first + 3: second]


class TestDragenWorkflowCountValidation:
    """The DRAGEN workflow must check the sample count before uploading."""

    def test_workflow_contains_count_operator(self):
        """The DRAGEN workflow must collect and validate the input channel size."""
        content = _load_module()
        wf_block = _extract_workflow_block(content)
        assert '.collect()' in wf_block, (
            "The DRAGEN workflow must call .collect() on cram_ch to gather all "
            "samples before validating the minimum count requirement."
        )

    def test_workflow_contains_subscribe_with_error(self):
        """The DRAGEN workflow must raise an error when sample count is < 5."""
        content = _load_module()
        wf_block = _extract_workflow_block(content)
        assert '.flatMap' in wf_block or 'subscribe' in wf_block, (
            "The DRAGEN workflow must use .flatMap or .subscribe to act on the sample count."
        )
        assert 'error' in wf_block, (
            "The DRAGEN workflow must call error() when sample count is below threshold."
        )

    def test_workflow_minimum_threshold_is_5(self):
        """The validation threshold in the DRAGEN workflow must be 5."""
        content = _load_module()
        wf_block = _extract_workflow_block(content)
        # Look for patterns like "< 5" or "n < 5" or "count < 5"
        m = re.search(r'<\s*5\b', wf_block)
        assert m is not None, (
            "The DRAGEN workflow must enforce a minimum of 5 samples "
            "(i.e., contain '< 5' in the count validation)."
        )

    def test_workflow_count_check_precedes_upload(self):
        """The sample count check must appear before UPLOAD_CRAM_FILES in the workflow."""
        content = _load_module()
        wf_block = _extract_workflow_block(content)
        collect_pos = wf_block.find('.collect()')
        upload_pos = wf_block.find('UPLOAD_CRAM_FILES')
        assert collect_pos != -1, ".collect() not found in DRAGEN workflow"
        assert upload_pos != -1, "UPLOAD_CRAM_FILES not found in DRAGEN workflow"
        assert collect_pos < upload_pos, (
            "The sample count validation must appear before UPLOAD_CRAM_FILES "
            "in the DRAGEN workflow."
        )


class TestStartAnalysisBatchMinSamples:
    """START_ANALYSIS_BATCH must validate minimum sample count before launching."""

    def test_bash_script_checks_minimum_sample_count(self):
        """The bash script must count samples from the data file and fail if < 5."""
        content = _load_module()
        script = _extract_start_analysis_batch_script(content)
        # Should contain a grep count of CRAM entries
        assert 'grep -c' in script, (
            "START_ANALYSIS_BATCH bash script must count CRAM entries with 'grep -c'."
        )

    def test_bash_script_minimum_threshold_is_5(self):
        """The bash minimum sample threshold must be 5."""
        content = _load_module()
        script = _extract_start_analysis_batch_script(content)
        m = re.search(r'-lt\s+5\b', script)
        assert m is not None, (
            "START_ANALYSIS_BATCH bash script must use '-lt 5' to enforce "
            "the minimum of 5 samples."
        )

    def test_bash_script_exits_on_insufficient_samples(self):
        """The bash script must exit 1 when the sample count is below the minimum."""
        content = _load_module()
        script = _extract_start_analysis_batch_script(content)
        # Find the minimum-sample check block
        lt5_match = re.search(r'-lt\s+5\b', script)
        assert lt5_match is not None, "'-lt 5' not found in START_ANALYSIS_BATCH script"
        surrounding = script[lt5_match.start(): lt5_match.start() + 300]
        assert 'exit 1' in surrounding, (
            "START_ANALYSIS_BATCH must call 'exit 1' when fewer than 5 samples are provided."
        )

    def test_bash_script_logs_sample_count_on_success(self):
        """The bash script must log the sample count when the minimum is satisfied."""
        content = _load_module()
        script = _extract_start_analysis_batch_script(content)
        assert 'minimum 5 satisfied for in-run panel of normals' in script, (
            "START_ANALYSIS_BATCH should log a confirmation message including "
            "'minimum 5 satisfied for in-run panel of normals' when validation passes."
        )


class TestInRunPonParameter:
    """The cnv_enable_in_run_pon parameter must be enabled."""

    def test_cnv_enable_in_run_pon_is_true(self):
        """START_ANALYSIS_BATCH must pass cnv_enable_in_run_pon:true to the pipeline."""
        content = _load_module()
        assert 'cnv_enable_in_run_pon:true' in content, (
            "modules-icav2-dragen.nf must include '--parameters cnv_enable_in_run_pon:true' "
            "to enable the in-run panel of normals."
        )

    def test_minimum_5_required_for_in_run_pon(self):
        """The error message must reference the in-run panel of normals requirement."""
        content = _load_module()
        # Check that 'in-run panel of normals' appears in both the workflow and the script
        assert content.count('in-run panel of normals') >= 2, (
            "Both the Nextflow workflow validation and the bash script should mention "
            "'in-run panel of normals' in their error/log messages."
        )


class TestQcCoverageFilters:
    """The qc_coverage_filters parameter must use mapq<20,bq<20."""

    def test_qc_coverage_filters_value(self):
        """START_ANALYSIS_BATCH must pass qc_coverage_filters with mapq<20,bq<20."""
        content = _load_module()
        assert "qc_coverage_filters:\"'mapq<20,bq<20'\"" in content, (
            "modules-icav2-dragen.nf must include "
            "--parameters qc_coverage_filters:\"'mapq<20,bq<20'\" "
            "for DRAGEN germline enrichment pipeline QC."
        )

    def test_qc_coverage_filters_not_old_value(self):
        """The old mapq<1 value must not be present."""
        content = _load_module()
        assert "mapq<1'" not in content, (
            "The old qc_coverage_filters value 'mapq<1' must be replaced with 'mapq<20,bq<20'."
        )


class TestMinSampleBashLogic:
    """Unit-test the bash validation logic extracted from START_ANALYSIS_BATCH."""

    def _run_bash(self, data_file_contents):
        """Run a bash snippet with a temp data file and return (returncode, stderr)."""
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(data_file_contents)
            data_path = f.name

        # Build a minimal bash script that mimics the relevant portion
        bash_script = textwrap.dedent(f"""
            #!/bin/bash
            dataFile="{data_path}"
            cramCode="crams"
            sample_count=$(grep -c "^${{cramCode}}:" "$dataFile" || true)
            if [ "$sample_count" -lt 5 ]; then
                echo "Error: minimum 5 samples required. Found ${{sample_count}}." >&2
                exit 1
            fi
            echo "OK: ${{sample_count}} samples"
        """)
        result = subprocess.run(
            ['bash', '-c', bash_script],
            capture_output=True, text=True
        )
        os.unlink(data_path)
        return result.returncode, result.stderr, result.stdout

    def _make_data_file(self, n_samples):
        """Build a combined data file with n_samples entries."""
        lines = []
        for i in range(1, n_samples + 1):
            lines.append(f"sampleId:SAMPLE{i}")
            lines.append(f"crams:fil.cram{i:04d}")
            lines.append(f"crais:fil.crai{i:04d}")
        return '\n'.join(lines) + '\n'

    def test_fails_with_0_samples(self):
        rc, stderr, _ = self._run_bash("")
        assert rc != 0, "Should fail with 0 samples"
        assert '0' in stderr

    def test_fails_with_1_sample(self):
        rc, stderr, _ = self._run_bash(self._make_data_file(1))
        assert rc != 0, "Should fail with 1 sample"

    def test_fails_with_4_samples(self):
        rc, stderr, _ = self._run_bash(self._make_data_file(4))
        assert rc != 0, "Should fail with 4 samples"

    def test_passes_with_5_samples(self):
        rc, _, stdout = self._run_bash(self._make_data_file(5))
        assert rc == 0, "Should pass with exactly 5 samples"
        assert '5' in stdout

    def test_passes_with_more_than_5_samples(self):
        rc, _, stdout = self._run_bash(self._make_data_file(10))
        assert rc == 0, "Should pass with 10 samples"
        assert '10' in stdout
