#!/usr/bin/env python3
"""
Tests for the normalise_cnv_caller_quality_scores.py script and its integration
into each CNV caller Nextflow module.

Validates:
  1. The script has a shebang, imports argparse and pysam, and exposes a CLI
     with --input_vcf, --output_vcf, and --caller arguments.
  2. Each caller is handled with the correct quality-score fields (FORMAT vs INFO).
  3. Each Nextflow module contains the NORMALISE_CNV_QUALITY_SCORES process with
     the correct --caller value.
  4. Each module's workflow calls the process and emits normalised_vcf /
     normalised_vcf_index.
  5. nextflow.config has a 'pysam' label entry for the container.
"""

import os
import re

import pytest

BIN_DIR  = os.path.join(os.path.dirname(__file__), '..', 'bin')
REPO_ROOT = os.path.join(os.path.dirname(__file__), '..')

SCRIPT_PATH = os.path.join(BIN_DIR, 'normalise_cnv_caller_quality_scores.py')

MODULES = {
    'CANOES':   'modules/callers/modules-canoes.nf',
    'CLAMMS':   'modules/callers/modules-clamms.nf',
    'XHMM':     'modules/callers/modules-xhmm.nf',
    'EXOMEDEPTH':'modules/callers/modules-exomedepth.nf',
    'GATK':     'modules/callers/modules-gatk-gcnv.nf',
    'CNVKIT':   'modules/callers/modules-cnvkit.nf',
    'DRAGEN':   'modules/callers/modules-icav2-dragen.nf',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_file(rel_path):
    path = os.path.join(REPO_ROOT, rel_path)
    with open(path) as fh:
        return fh.read()


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _extract_process(module_text, process_name):
    match = re.search(
        rf"process {re.escape(process_name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def _uses_shared_process(module_text, process_name='NORMALISE_CNV_QUALITY_SCORES'):
    """Return True when the module imports process_name from modules-common.nf."""
    return bool(re.search(
        rf"include\s*\{{[^}}]*{re.escape(process_name)}[^}}]*\}}.*modules-common\.nf",
        module_text,
    ))


def _extract_workflow(module_text, workflow_name):
    match = re.search(
        rf"workflow {re.escape(workflow_name)}\s*\{{(.+)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def canoes_text():
    return _read_file('modules/callers/modules-canoes.nf')


@pytest.fixture(scope='module')
def clamms_text():
    return _read_file('modules/callers/modules-clamms.nf')


@pytest.fixture(scope='module')
def xhmm_text():
    return _read_file('modules/callers/modules-xhmm.nf')


@pytest.fixture(scope='module')
def exomedepth_text():
    return _read_file('modules/callers/modules-exomedepth.nf')


@pytest.fixture(scope='module')
def gatk_text():
    return _read_file('modules/callers/modules-gatk-gcnv.nf')


@pytest.fixture(scope='module')
def cnvkit_text():
    return _read_file('modules/callers/modules-cnvkit.nf')


@pytest.fixture(scope='module')
def dragen_text():
    return _read_file('modules/callers/modules-icav2-dragen.nf')


@pytest.fixture(scope='module')
def config_text():
    return _read_file('nextflow.config')


# ===========================================================================
# 1. Script structure and CLI
# ===========================================================================

class TestScriptStructure:
    """normalise_cnv_caller_quality_scores.py must be a well-formed CLI script."""

    def test_shebang_present(self, script_text):
        assert script_text.startswith('#!/usr/bin/env python3'), (
            "Script must start with #!/usr/bin/env python3 shebang"
        )

    def test_imports_argparse(self, script_text):
        assert 'import argparse' in script_text, (
            "Script must import argparse to provide CLI argument parsing"
        )

    def test_imports_pysam(self, script_text):
        assert 'import pysam' in script_text, (
            "Script must import pysam to read and write VCF files"
        )

    def test_has_main_block(self, script_text):
        assert "__name__ == '__main__'" in script_text or \
               '__name__ == "__main__"' in script_text, (
            "Script must have an if __name__ == '__main__' guard so it is "
            "usable both as a module and as a standalone CLI tool"
        )

    def test_cli_has_input_vcf(self, script_text):
        assert '--input_vcf' in script_text, (
            "CLI must accept --input_vcf argument"
        )

    def test_cli_has_output_vcf(self, script_text):
        assert '--output_vcf' in script_text, (
            "CLI must accept --output_vcf argument"
        )

    def test_cli_has_caller(self, script_text):
        assert '--caller' in script_text, (
            "CLI must accept --caller argument to identify the CNV caller"
        )

    def test_cli_choices_include_all_callers(self, script_text):
        for caller in ('CANOES', 'CLAMMS', 'XHMM', 'GATK', 'DRAGEN', 'CNVKIT', 'EXOMEDEPTH'):
            assert caller in script_text, (
                f"CLI --caller choices must include {caller}"
            )


# ===========================================================================
# 2. Caller-specific quality-score field correctness
# ===========================================================================

class TestCallerFieldCorrectness:
    """Each caller block must read quality scores from the correct VCF location."""

    def test_canoes_reads_from_format_not_info(self, script_text):
        """CANOES stores Q_SOME in FORMAT, not INFO."""
        canoes_block = re.search(
            r'if caller == "CANOES":(.+?)(?=elif caller|$)',
            script_text, re.DOTALL
        )
        assert canoes_block, "CANOES caller block not found"
        block = canoes_block.group(1)
        # Must use sample["Q_SOME"] (FORMAT), not record.info["Q_SOME"]
        assert 'record.info' not in block, (
            "CANOES block must not access Q_SOME via record.info; "
            "Q_SOME is a FORMAT field in CANOES VCFs"
        )
        assert 'sample' in block and '"Q_SOME"' in block, (
            "CANOES block must read Q_SOME from sample FORMAT fields"
        )

    def test_clamms_reads_from_format_not_info(self, script_text):
        """CLAMMS stores Q_SOME and Q_EXACT in FORMAT, not INFO."""
        clamms_block = re.search(
            r'elif caller == "CLAMMS":(.+?)(?=elif caller|$)',
            script_text, re.DOTALL
        )
        assert clamms_block, "CLAMMS caller block not found"
        block = clamms_block.group(1)
        # Must use sample["Q_SOME"] and sample["Q_EXACT"] (FORMAT), not record.info
        assert 'record.info' not in block, (
            "CLAMMS block must not access Q_SOME/Q_EXACT via record.info; "
            "they are FORMAT fields in CLAMMS VCFs"
        )
        assert 'sample' in block and '"Q_SOME"' in block, (
            "CLAMMS block must read Q_SOME from sample FORMAT fields"
        )
        assert '"Q_EXACT"' in block, (
            "CLAMMS block must read Q_EXACT from sample FORMAT fields "
            "to gate the normalisation: Q_EXACT < 0 → QUAL_norm = 0"
        )

    def test_sample_accessor_uses_first_value(self, script_text):
        """sample must be obtained from record.samples.values(), not record.samples."""
        assert 'record.samples.values()' in script_text, (
            "sample must be obtained via next(iter(record.samples.values())) "
            "to correctly access the single-sample FORMAT fields in pysam"
        )
        # The old broken assignment was: sample = record.samples
        assert 'sample = record.samples\n' not in script_text and \
               'sample = record.samples  ' not in script_text, (
            "Assigning sample = record.samples is incorrect; it returns the "
            "VariantRecordSamples dict-like object, not a single sample's FORMAT fields"
        )


# ===========================================================================
# 3. NORMALISE_CNV_QUALITY_SCORES process defined in modules-common.nf and wired in main.nf
# ===========================================================================

class TestNormaliseProcessExists:
    """NORMALISE_CNV_QUALITY_SCORES must be defined in modules-common.nf and
    modules-normalise.nf, and imported into main.nf for each caller."""

    @pytest.fixture(scope='class')
    def common_text(self):
        return _read_file('modules/common/modules-common.nf')

    @pytest.fixture(scope='class')
    def main_text(self):
        return _read_file('main.nf')

    def test_common_defines_normalise_process(self, common_text):
        body = _extract_process(common_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert body is not None, "NORMALISE_CNV_QUALITY_SCORES process not found in modules-common.nf"

    def test_main_imports_normalise_for_all_callers(self, main_text):
        for caller in ['CANOES', 'CLAMMS', 'XHMM', 'CNVKIT', 'GCNV', 'DRAGEN', 'EXOMEDEPTH']:
            assert f"NORM_{caller}_FULL" in main_text, (
                f"main.nf must import or alias NORMALISE_CNV_QUALITY_SCORES as NORM_{caller}_FULL"
            )


class TestNormaliseProcessCaller:
    """main.nf must pass the correct --caller value to NORM_<CALLER>_FULL after filtering."""

    @pytest.fixture(scope='class')
    def main_text(self):
        return _read_file('main.nf')

    def _check_caller(self, main_text, norm_name, expected_caller):
        pattern = rf"{norm_name}\([^)]*'{expected_caller}'"
        assert re.search(pattern, main_text), (
            f"main.nf must invoke {norm_name} passing '{expected_caller}' as the caller argument"
        )

    def test_canoes_caller(self, main_text):
        self._check_caller(main_text, 'NORM_CANOES_FULL', 'CANOES')

    def test_clamms_caller(self, main_text):
        self._check_caller(main_text, 'NORM_CLAMMS_FULL', 'CLAMMS')

    def test_xhmm_caller(self, main_text):
        self._check_caller(main_text, 'NORM_XHMM_FULL', 'XHMM')

    def test_exomedepth_caller(self, main_text):
        self._check_caller(main_text, 'NORM_EXOMEDEPTH_FULL', 'EXOMEDEPTH')

    def test_gatk_caller(self, main_text):
        self._check_caller(main_text, 'NORM_GCNV_FULL', 'GATK')

    def test_cnvkit_caller(self, main_text):
        self._check_caller(main_text, 'NORM_CNVKIT_FULL', 'CNVKIT')

    def test_dragen_caller(self, main_text):
        self._check_caller(main_text, 'NORM_DRAGEN_FULL', 'DRAGEN')


class TestNormaliseProcessOutputs:
    """NORMALISE_CNV_QUALITY_SCORES must emit normalised_vcf and normalised_vcf_index."""

    @pytest.fixture(scope='class')
    def common_text(self):
        return _read_file('modules/common/modules-common.nf')

    def test_common_process_emits_normalised_vcf(self, common_text):
        body = _extract_process(common_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert body is not None, "NORMALISE_CNV_QUALITY_SCORES process not found in modules-common.nf"
        assert 'normalised_vcf' in body, (
            "NORMALISE_CNV_QUALITY_SCORES must emit normalised_vcf"
        )
        assert 'normalised_vcf_index' in body, (
            "NORMALISE_CNV_QUALITY_SCORES must emit normalised_vcf_index"
        )
        assert '.normalised' in body, (
            "Output files must use .normalised VCF suffix"
        )


# ===========================================================================
# 4. Caller workflows emit raw calls; normalisation is deferred until after filtering
# ===========================================================================

class TestCallerWorkflowsDoNotPreNormalise:
    """Caller modules must not normalise raw calls prematurely.
    
    Quality score normalisation must strictly occur after quality score filtering
    (in main.nf and standalone post-filtering workflows).
    """

    def _check_workflow_does_not_call_normalise(self, module_text, workflow_name, caller_name):
        wf = _extract_workflow(module_text, workflow_name)
        assert wf is not None, f"{workflow_name} workflow not found"
        assert 'NORMALISE_CNV_QUALITY_SCORES' not in wf, (
            f"{caller_name} caller workflow must NOT call NORMALISE_CNV_QUALITY_SCORES on raw calls; "
            "normalisation must only happen after quality score filtering"
        )

    def test_canoes_no_prenormalise(self, canoes_text):
        self._check_workflow_does_not_call_normalise(canoes_text, 'CANOES', 'CANOES')

    def test_clamms_no_prenormalise(self, clamms_text):
        self._check_workflow_does_not_call_normalise(clamms_text, 'CLAMMS', 'CLAMMS')

    def test_xhmm_no_prenormalise(self, xhmm_text):
        self._check_workflow_does_not_call_normalise(xhmm_text, 'XHMM', 'XHMM')

    def test_exomedepth_no_prenormalise(self, exomedepth_text):
        self._check_workflow_does_not_call_normalise(exomedepth_text, 'EXOMEDEPTH', 'EXOMEDEPTH')

    def test_gatk_no_prenormalise(self, gatk_text):
        self._check_workflow_does_not_call_normalise(gatk_text, 'GATK_GCNV', 'GATK')

    def test_cnvkit_no_prenormalise(self, cnvkit_text):
        self._check_workflow_does_not_call_normalise(cnvkit_text, 'CNVKIT', 'CNVKIT')

    def test_dragen_no_prenormalise(self, dragen_text):
        self._check_workflow_does_not_call_normalise(dragen_text, 'DRAGEN', 'DRAGEN')


class TestWorkflowEmitsRawSorted:
    """Each caller workflow emit block must expose raw sorted_vcf and sorted_vcf_index."""

    def _check_emit(self, module_text, workflow_name, caller_name):
        wf = _extract_workflow(module_text, workflow_name)
        assert wf is not None, f"{workflow_name} workflow not found"
        assert 'sorted_vcf' in wf, (
            f"{caller_name} workflow must emit sorted_vcf"
        )
        assert 'sorted_vcf_index' in wf, (
            f"{caller_name} workflow must emit sorted_vcf_index"
        )

    def test_canoes_emit(self, canoes_text):
        self._check_emit(canoes_text, 'CANOES', 'CANOES')

    def test_clamms_emit(self, clamms_text):
        self._check_emit(clamms_text, 'CLAMMS', 'CLAMMS')

    def test_xhmm_emit(self, xhmm_text):
        self._check_emit(xhmm_text, 'XHMM', 'XHMM')

    def test_exomedepth_emit(self, exomedepth_text):
        self._check_emit(exomedepth_text, 'EXOMEDEPTH', 'EXOMEDEPTH')

    def test_gatk_emit(self, gatk_text):
        self._check_emit(gatk_text, 'GATK_GCNV', 'GATK')

    def test_cnvkit_emit(self, cnvkit_text):
        self._check_emit(cnvkit_text, 'CNVKIT', 'CNVKIT')

    def test_dragen_emit(self, dragen_text):
        self._check_emit(dragen_text, 'DRAGEN', 'DRAGEN')


# ===========================================================================
# 5. nextflow.config has pysam container label
# ===========================================================================

class TestNextflowConfigPysamLabel:
    """nextflow.config must declare a 'pysam' label for the container."""

    def test_pysam_label_present(self, config_text):
        assert "withLabel: 'pysam'" in config_text, (
            "nextflow.config must contain a withLabel: 'pysam' block "
            "so that NORMALISE_CNV_QUALITY_SCORES processes can resolve a container"
        )

    def test_pysam_container_uses_pysam_image(self, config_text):
        pysam_block_match = re.search(
            r"withLabel:\s*'pysam'\s*\{(.+?)(?=\})", config_text, re.DOTALL
        )
        assert pysam_block_match, "pysam label block not found in nextflow.config"
        block = pysam_block_match.group(1)
        assert 'pysam' in block.lower(), (
            "The 'pysam' container image path/URI must reference pysam"
        )

    def test_pysam_container_uses_expected_biocontainers_tag(self, config_text):
        pysam_block_match = re.search(
            r"withLabel:\s*'pysam'\s*\{(.+?)(?=\})", config_text, re.DOTALL
        )
        assert pysam_block_match, "pysam label block not found in nextflow.config"
        block = pysam_block_match.group(1)
        assert (
            "docker://quay.io/biocontainers/pysam:0.23.3--py313hd07c5dd_2" in block
        ), "The pysam label must use the expected biocontainers pysam tag"


# ===========================================================================
# 6. DRAGEN QUAL normalisation chain
# ===========================================================================

class TestDragenQualNormalisationChain:
    """Verify the DRAGEN-specific QUAL normalisation chain in the script.

    DRAGEN CNV VCFs produced by DRAGEN Germline Enrichment carry native QUAL
    values (QUAL_{native}) in the QUAL column and have not been INFO-tagged.
    After the processing chain the VCF must:
      1. Be INFO-tagged with TOOL=DRAGEN.
      2. Be bgzipped, bcftools-sorted, and tabix-indexed.
      3. Have QUAL replaced with the normalised score (Q_norm).
      4. Have the original QUAL preserved in the FORMAT/SAMPLE field OQ.
      5. Have OAS set to "QUAL" to record which metric was normalised.
    """

    def _dragen_block(self, script_text):
        m = re.search(
            r'elif caller == "DRAGEN":(.+?)(?=\nelif caller |\nexcept |\Z)',
            script_text, re.DOTALL,
        )
        assert m, "DRAGEN caller block not found in normalise_cnv_caller_quality_scores.py"
        return m.group(1)

    # -------------------------------------------------------------------------
    # OQ / OAS FORMAT headers must be declared so pysam can write them
    # -------------------------------------------------------------------------

    def test_oq_format_header_declared(self, script_text):
        """Script must declare an OQ FORMAT header to store the original QUAL."""
        assert '"OQ"' in script_text or "'OQ'" in script_text, (
            "normalise_cnv_caller_quality_scores.py must declare "
            "vcf.header.formats.add('OQ', ...) to carry the original QUAL value "
            "in the FORMAT/SAMPLE field of every normalised record"
        )

    def test_oas_format_header_declared(self, script_text):
        """Script must declare an OAS FORMAT header to store the metric name."""
        assert '"OAS"' in script_text or "'OAS'" in script_text, (
            "normalise_cnv_caller_quality_scores.py must declare "
            "vcf.header.formats.add('OAS', ...) to record which metric was used"
        )

    # -------------------------------------------------------------------------
    # Original QUAL is moved to OQ before being cleared
    # -------------------------------------------------------------------------

    def test_original_qual_moved_to_oq_format_field(self, script_text):
        """The script must store orig_qual in sample['OQ'] before clearing QUAL."""
        assert 'sample["OQ"]' in script_text or "sample['OQ']" in script_text, (
            "normalise_cnv_caller_quality_scores.py must assign sample['OQ'] = orig_qual "
            "so the original QUAL_{native} is preserved in the FORMAT/SAMPLE field"
        )

    def test_qual_cleared_before_normalisation(self, script_text):
        """QUAL must be explicitly cleared (set to None) before writing the normalised value."""
        assert 'record.qual = None' in script_text, (
            "normalise_cnv_caller_quality_scores.py must set record.qual = None "
            "after saving orig_qual to OQ, so the QUAL field is properly replaced "
            "with the normalised score"
        )

    def test_normalised_qual_replaces_qual_field(self, script_text):
        """record.qual must be assigned the normalised score after the caller block."""
        assert 'record.qual = round(qual_norm' in script_text, (
            "normalise_cnv_caller_quality_scores.py must assign "
            "record.qual = round(qual_norm, ...) to replace the native QUAL with Q_norm"
        )

    # -------------------------------------------------------------------------
    # DRAGEN block reads from record.qual, not FORMAT fields
    # -------------------------------------------------------------------------

    def test_dragen_reads_from_record_qual_not_format(self, script_text):
        """DRAGEN block must use orig_qual (from record.qual), not sample FORMAT fields."""
        block = self._dragen_block(script_text)
        # orig_qual was captured from record.qual before the caller section
        assert 'orig_qual' in block, (
            "DRAGEN caller block must read from orig_qual (i.e. record.qual), "
            "not from sample FORMAT fields; DRAGEN stores the quality score in QUAL"
        )

    def test_dragen_does_not_read_format_fields_for_qual(self, script_text):
        """DRAGEN block must not extract quality from a FORMAT field."""
        block = self._dragen_block(script_text)
        # DRAGEN QUAL is in the record QUAL column, not in FORMAT subfields like Q_SOME
        for field in ('"Q_SOME"', '"SQ"', '"CNQ"', '"QS"'):
            assert field not in block, (
                f"DRAGEN caller block must not use FORMAT field {field} for quality "
                "normalisation; DRAGEN stores the score in the QUAL column"
            )

    # -------------------------------------------------------------------------
    # DRAGEN metric_used = "QUAL" → written to OAS
    # -------------------------------------------------------------------------

    def test_dragen_metric_used_is_qual(self, script_text):
        """DRAGEN block must set metric_used = 'QUAL'."""
        block = self._dragen_block(script_text)
        assert 'metric_used = "QUAL"' in block or "metric_used = 'QUAL'" in block, (
            "DRAGEN caller block must set metric_used = 'QUAL' so that OAS records "
            "which score was used for normalisation"
        )

    # -------------------------------------------------------------------------
    # DRAGEN only normalises PASS records
    # -------------------------------------------------------------------------

    def test_dragen_gates_on_pass_filter(self, script_text):
        """DRAGEN block must only normalise records whose FILTER contains PASS."""
        block = self._dragen_block(script_text)
        assert 'PASS' in block, (
            "DRAGEN caller block must check that the record has a PASS filter "
            "before normalising QUAL; non-PASS records receive qual_norm = 0.0"
        )
        assert 'record.filter' in block, (
            "DRAGEN caller block must inspect record.filter to determine "
            "whether the record has passed variant-calling filters"
        )
