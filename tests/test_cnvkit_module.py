#!/usr/bin/env python3
"""
Tests for the CNVKit Nextflow module (modules/callers/modules-cnvkit.nf).

Validates that the module correctly implements the CNVKit pipeline for exome
sequencing data with a cohort-only sampleset (no matched normals), producing
per-sample VCFs:

  1. CALL_CNV runs `cnvkit.py call` after `cnvkit.py segment` to produce
     integer copy-number calls (.call.cns). Without this step, `export vcf`
     lacks the CN column and cannot produce a proper CNV VCF.

  2. CALL_CNV emits the called-segments file (.call.cns) in its output tuple
     so that EXPORT_RESULTS can consume it.

  3. EXPORT_RESULTS accepts a 4-element tuple (sample_id, cnr, cns, call_cns)
     and uses the called-segments file ($call_cns) -- not the raw segmented
     file ($cns) -- for both `export vcf` and `export bed`.

  4. The pooled-reference approach (BUILD from all cohort coverage files) is
     intact, confirming cohort-mode operation without matched normals.
"""

import os
import re

import pytest


NF_MODULE = "modules/callers/modules-cnvkit.nf"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_module(repo_root="."):
    path = os.path.join(repo_root, NF_MODULE)
    with open(path) as fh:
        return fh.read()


def _extract_process(module_text, process_name):
    """Return the body of a named process block, or None if not found."""
    match = re.search(
        rf"process {re.escape(process_name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def module_text():
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    return _read_module(repo_root)


class TestParamDeclarationStrategy:
    """CNVkit module should consume globally-declared params."""

    def test_module_has_no_strict_params_block(self, module_text):
        assert "params {" not in module_text, (
            "modules-cnvkit.nf should not define a strict module-level params block; "
            "pipeline-facing params are centralized globally."
        )


@pytest.fixture(scope="module")
def call_cnv_body(module_text):
    body = _extract_process(module_text, "CALL_CNV")
    assert body is not None, "CALL_CNV process not found in module"
    return body


@pytest.fixture(scope="module")
def export_results_body(module_text):
    body = _extract_process(module_text, "EXPORT_RESULTS")
    assert body is not None, "EXPORT_RESULTS process not found in module"
    return body


# ===========================================================================
# 1. CALL_CNV: cnvkit.py call must be present after cnvkit.py segment
# ===========================================================================

class TestCallCnvProcess:
    """CALL_CNV must run `cnvkit.py call` to produce discrete integer CN calls."""

    def test_segment_command_present(self, call_cnv_body):
        """cnvkit.py segment must be present (produces log2-ratio segments)."""
        assert "cnvkit.py segment" in call_cnv_body, (
            "cnvkit.py segment not found in CALL_CNV process"
        )

    def test_call_command_present(self, call_cnv_body):
        """cnvkit.py call must follow segment to assign integer copy numbers.

        Without this step the .cns file has no 'cn' column and
        `cnvkit.py export vcf` cannot produce a valid CNV VCF.
        """
        assert "cnvkit.py call" in call_cnv_body, (
            "cnvkit.py call is missing from CALL_CNV. "
            "It must run after cnvkit.py segment to produce integer CN calls "
            "(.call.cns) required by cnvkit.py export vcf."
        )

    def test_call_precedes_export_in_script(self, call_cnv_body):
        """cnvkit.py call must appear after cnvkit.py segment in the script."""
        seg_pos = call_cnv_body.find("cnvkit.py segment")
        call_pos = call_cnv_body.find("cnvkit.py call")
        assert seg_pos != -1 and call_pos != -1, (
            "Both cnvkit.py segment and cnvkit.py call must be present"
        )
        assert call_pos > seg_pos, (
            "cnvkit.py call must appear after cnvkit.py segment in the script"
        )

    def test_output_emits_call_cns(self, call_cnv_body):
        """CALL_CNV output tuple must include the .call.cns file."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", call_cnv_body, re.DOTALL
        )
        assert output_section, "output section not found in CALL_CNV"
        out_text = output_section.group(1)
        assert ".call.cns" in out_text, (
            "CALL_CNV output must emit the .call.cns file (from cnvkit.py call) "
            "so that EXPORT_RESULTS can use it for VCF/BED export."
        )

    def test_output_is_four_element_tuple(self, call_cnv_body):
        """CALL_CNV results tuple must contain 4 elements: sample_id, cnr, cns, call_cns."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", call_cnv_body, re.DOTALL
        )
        assert output_section, "output section not found in CALL_CNV"
        out_text = output_section.group(1)
        # Expect a tuple with val(sample_id), path(.cnr), path(.cns), path(.call.cns)
        assert re.search(r"tuple\s+val\(sample_id\)", out_text), (
            "CALL_CNV results emit must start with val(sample_id)"
        )
        assert ".cnr" in out_text, "CALL_CNV results emit must include .cnr"
        assert ".cns" in out_text, "CALL_CNV results emit must include .cns"
        assert ".call.cns" in out_text, "CALL_CNV results emit must include .call.cns"


# ===========================================================================
# 2. EXPORT_RESULTS: must use called segments (.call.cns) for VCF/BED export
# ===========================================================================

class TestExportResultsProcess:
    """EXPORT_RESULTS must consume the called-segments file for VCF/BED output."""

    def test_input_accepts_call_cns(self, export_results_body):
        """Input tuple must include a path for the called-segments file."""
        input_section = re.search(
            r"input:\s*(.+?)(?=\s*output:)", export_results_body, re.DOTALL
        )
        assert input_section, "input section not found in EXPORT_RESULTS"
        in_text = input_section.group(1)
        assert "call_cns" in in_text, (
            "EXPORT_RESULTS input tuple must include call_cns (the .call.cns file "
            "produced by cnvkit.py call) so that integer CN values are available "
            "for cnvkit.py export vcf."
        )

    def test_input_is_four_element_tuple(self, export_results_body):
        """Input tuple must be (sample_id, cnr, cns, call_cns)."""
        input_section = re.search(
            r"input:\s*(.+?)(?=\s*output:)", export_results_body, re.DOTALL
        )
        assert input_section, "input section not found in EXPORT_RESULTS"
        in_text = input_section.group(1)
        assert re.search(r"tuple\s+val\(sample_id\)", in_text), (
            "EXPORT_RESULTS input tuple must start with val(sample_id)"
        )
        assert "cnr" in in_text, "EXPORT_RESULTS input tuple must include cnr"
        assert "cns" in in_text, "EXPORT_RESULTS input tuple must include cns"
        assert "call_cns" in in_text, "EXPORT_RESULTS input tuple must include call_cns"

    def test_export_vcf_uses_call_cns(self, export_results_body):
        """cnvkit.py export vcf must use $call_cns (not $cns).

        $cns (from segment) lacks integer CN values; only $call_cns (from call)
        has the 'cn' column required for a valid VCF export.
        """
        script_section = re.search(
            r'"""(.+?)"""', export_results_body, re.DOTALL
        )
        assert script_section, "script block not found in EXPORT_RESULTS"
        script_text = script_section.group(1)

        # Find the export vcf line
        vcf_line_match = re.search(r"cnvkit\.py export vcf(.+)", script_text)
        assert vcf_line_match, "cnvkit.py export vcf not found in EXPORT_RESULTS script"
        vcf_line = vcf_line_match.group(1)

        assert "$call_cns" in vcf_line, (
            "cnvkit.py export vcf must use $call_cns (the output of cnvkit.py call "
            "with integer CN values), not $cns (the raw segmented output)."
        )

    def test_export_bed_uses_call_cns(self, export_results_body):
        """cnvkit.py export bed must use $call_cns (not $cns)."""
        script_section = re.search(
            r'"""(.+?)"""', export_results_body, re.DOTALL
        )
        assert script_section, "script block not found in EXPORT_RESULTS"
        script_text = script_section.group(1)

        bed_line_match = re.search(r"cnvkit\.py export bed(.+)", script_text)
        assert bed_line_match, "cnvkit.py export bed not found in EXPORT_RESULTS script"
        bed_line = bed_line_match.group(1)

        assert "$call_cns" in bed_line, (
            "cnvkit.py export bed must use $call_cns (called segments with integer CN), "
            "not $cns (raw segmented output)."
        )


# ===========================================================================
# 3. Cohort-mode: pooled reference built from all sample coverages
# ===========================================================================

class TestPooledReferenceApproach:
    """CREATE_POOLED_REFERENCE must build a reference from ALL cohort coverages."""

    def test_pooled_reference_process_exists(self, module_text):
        """CREATE_POOLED_REFERENCE process must be present for cohort-mode operation."""
        assert "process CREATE_POOLED_REFERENCE" in module_text, (
            "CREATE_POOLED_REFERENCE process not found. "
            "A pooled reference built from all cohort samples is required when "
            "running CNVKit on a cohort without matched normals."
        )

    def test_reference_command_uses_all_cnn(self, module_text):
        """cnvkit.py reference must consume coverage files from all cohort samples."""
        ref_body = _extract_process(module_text, "CREATE_POOLED_REFERENCE")
        assert ref_body is not None, "CREATE_POOLED_REFERENCE process not found"
        script_section = re.search(r'"""(.+?)"""', ref_body, re.DOTALL)
        assert script_section, "script block not found in CREATE_POOLED_REFERENCE"
        script_text = script_section.group(1)
        assert "cnvkit.py reference" in script_text, (
            "cnvkit.py reference must be used to build a pooled reference"
        )
        assert "*.cnn" in script_text, (
            "cnvkit.py reference must consume all .cnn coverage files (*.cnn) "
            "from the cohort to build an unbiased pooled reference."
        )

    def test_workflow_collects_all_coverages(self, module_text):
        """CNVKIT workflow must collect coverage from all samples before building reference."""
        wf_match = re.search(
            r"workflow CNVKIT\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "CNVKIT workflow not found"
        wf_body = wf_match.group(1)
        assert ".collect()" in wf_body, (
            "CNVKIT workflow must use .collect() to gather all sample coverages "
            "before passing them to CREATE_POOLED_REFERENCE."
        )
