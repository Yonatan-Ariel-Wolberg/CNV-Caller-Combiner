#!/usr/bin/env python3
"""
Tests for DRAGEN TOOL annotation consistency and BGZIP_SORT_INDEX_VCF.

Validates that the DRAGEN module:
  1. Uses INFO field key TOOL (singular), consistent with all other 6 callers.
  2. Has a BGZIP_SORT_INDEX_VCF process to produce sorted, tabix-indexed VCFs.
  3. Chains annotation → sort/index → normalise in the DRAGEN workflow.
  4. Emits sorted_vcf and sorted_vcf_index from the DRAGEN workflow.
5. main.nf RUN_DRAGEN uses DRAGEN.out.normalised_vcf for EVALUATE so
   downstream checks use QUAL_norm from QUAL (consistent with other callers).
"""

import os
import re

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DRAGEN_MODULE = os.path.join(REPO_ROOT, 'modules', 'callers', 'modules-icav2-dragen.nf')
MAIN_NF = os.path.join(REPO_ROOT, 'main.nf')


def _load(path):
    with open(path) as f:
        return f.read()


def _extract_process(text, name):
    m = re.search(
        rf'process {re.escape(name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)',
        text, re.DOTALL)
    return m.group(1) if m else None


def _extract_workflow(text, name):
    m = re.search(
        rf'workflow {re.escape(name)}\s*\{{(.+)',
        text, re.DOTALL)
    return m.group(1) if m else None


@pytest.fixture(scope='module')
def dragen_text():
    return _load(DRAGEN_MODULE)


@pytest.fixture(scope='module')
def main_text():
    return _load(MAIN_NF)


# ---------------------------------------------------------------------------
# 1. ADD_DRAGEN_TOOL_ANNOTATION uses INFO/TOOL (singular)
# ---------------------------------------------------------------------------

class TestDragenToolAnnotationInfoField:
    """ADD_DRAGEN_TOOL_ANNOTATION must use TOOL (singular) in the INFO field."""

    def test_info_id_is_tool_singular(self, dragen_text):
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found in modules-icav2-dragen.nf"
        assert 'ID=TOOL,' in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must declare ##INFO=<ID=TOOL,...> "
            "(singular), consistent with all other callers"
        )

    def test_no_info_id_tools_plural(self, dragen_text):
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found"
        assert 'ID=TOOLS,' not in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must NOT use ##INFO=<ID=TOOLS,...> (plural); "
            "use TOOL (singular) to be consistent with all other callers"
        )

    def test_bcftools_annotate_uses_info_tool(self, dragen_text):
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found"
        assert 'INFO/TOOL' in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must use -c CHROM,FROM,TO,INFO/TOOL "
            "(singular) in bcftools annotate command"
        )
        assert 'INFO/TOOLS' not in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must NOT use INFO/TOOLS (plural)"
        )

    def test_annotated_output_uses_dragen_suffix_not_dragen_output(self, dragen_text):
        """Output filename must use _DRAGEN.annotated.vcf.gz (not _DRAGEN_output).

        The convention for all callers is ${sample_id}_CALLERNAME.normalised.vcf.gz.
        Using _DRAGEN (not _DRAGEN_output) as the intermediate suffix ensures the
        final normalised VCF is named ${sample_id}_DRAGEN.normalised.vcf.gz,
        matching what main.nf build_feature_inputs expects.
        """
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found"
        assert '_DRAGEN.annotated.vcf.gz' in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION output must use *_DRAGEN.annotated.vcf.gz "
            "so the downstream chain produces ${sample_id}_DRAGEN.normalised.vcf.gz"
        )
        assert '_DRAGEN_output.annotated.vcf.gz' not in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must NOT use *_DRAGEN_output.annotated.vcf.gz; "
            "this would produce _DRAGEN_output.normalised.vcf.gz instead of "
            "_DRAGEN.normalised.vcf.gz which main.nf build_feature_inputs expects"
        )


# ---------------------------------------------------------------------------
# 2. BGZIP_SORT_INDEX_VCF process exists in DRAGEN module
# ---------------------------------------------------------------------------

class TestDragenBgzipSortIndexVcf:
    """DRAGEN module must contain a BGZIP_SORT_INDEX_VCF process."""

    def test_process_exists(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, (
            "BGZIP_SORT_INDEX_VCF process not found in modules-icav2-dragen.nf; "
            "DRAGEN must sort and index VCFs for consistency with other callers"
        )

    def test_emits_sorted_vcf(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert 'sorted_vcf' in process_body, (
            "BGZIP_SORT_INDEX_VCF must emit sorted_vcf"
        )

    def test_emits_sorted_vcf_index(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert 'sorted_vcf_index' in process_body, (
            "BGZIP_SORT_INDEX_VCF must emit sorted_vcf_index"
        )

    def test_uses_sorted_vcf_gz_suffix(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert '.sorted.vcf.gz' in process_body, (
            "BGZIP_SORT_INDEX_VCF must produce *.sorted.vcf.gz files"
        )

    def test_uses_bcftools_sort(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert 'bcftools sort' in process_body, (
            "BGZIP_SORT_INDEX_VCF must use bcftools sort to sort the VCF"
        )

    def test_uses_tabix_indexing(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert 'tabix' in process_body, (
            "BGZIP_SORT_INDEX_VCF must use tabix to index the sorted VCF"
        )

    def test_label_is_bcftools(self, dragen_text):
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert "label 'bcftools'" in process_body, (
            "BGZIP_SORT_INDEX_VCF must use label 'bcftools' (consistent with "
            "the same process in all other caller modules)"
        )


# ---------------------------------------------------------------------------
# 3. DRAGEN workflow chains annotation → sort/index → normalise
# ---------------------------------------------------------------------------

class TestDragenWorkflowChaining:
    """DRAGEN workflow must chain BGZIP_SORT_INDEX_VCF between annotation and normalisation."""

    def test_bgzip_called_in_workflow(self, dragen_text):
        wf = _extract_workflow(dragen_text, 'DRAGEN')
        assert wf is not None, "DRAGEN workflow not found"
        assert 'BGZIP_SORT_INDEX_VCF' in wf, (
            "DRAGEN workflow must call BGZIP_SORT_INDEX_VCF"
        )

    def test_normalise_uses_sorted_vcf(self, dragen_text):
        wf = _extract_workflow(dragen_text, 'DRAGEN')
        assert wf is not None, "DRAGEN workflow not found"
        assert 'BGZIP_SORT_INDEX_VCF.out.sorted_vcf' in wf, (
            "DRAGEN workflow must pass BGZIP_SORT_INDEX_VCF.out.sorted_vcf "
            "to NORMALISE_CNV_QUALITY_SCORES"
        )

    def test_workflow_emits_sorted_vcf(self, dragen_text):
        wf = _extract_workflow(dragen_text, 'DRAGEN')
        assert wf is not None, "DRAGEN workflow not found"
        assert 'sorted_vcf' in wf, (
            "DRAGEN workflow emit block must include sorted_vcf"
        )

    def test_workflow_emits_sorted_vcf_index(self, dragen_text):
        wf = _extract_workflow(dragen_text, 'DRAGEN')
        assert wf is not None, "DRAGEN workflow not found"
        assert 'sorted_vcf_index' in wf, (
            "DRAGEN workflow emit block must include sorted_vcf_index"
        )


# ---------------------------------------------------------------------------
# 4. NORMALISE_CNV_QUALITY_SCORES strips .sorted.vcf.gz
# ---------------------------------------------------------------------------

class TestDragenNormaliseInputSuffix:
    """DRAGEN NORMALISE_CNV_QUALITY_SCORES must strip .sorted.vcf.gz (not .annotated.vcf.gz)."""

    def test_strips_sorted_suffix(self, dragen_text):
        process_body = _extract_process(dragen_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert process_body is not None, \
            "NORMALISE_CNV_QUALITY_SCORES process not found in modules-icav2-dragen.nf"
        assert "'.sorted.vcf.gz'" in process_body or '".sorted.vcf.gz"' in process_body, (
            "NORMALISE_CNV_QUALITY_SCORES in DRAGEN module must strip .sorted.vcf.gz "
            "from the input filename to derive the sample name"
        )

    def test_does_not_strip_annotated_suffix(self, dragen_text):
        process_body = _extract_process(dragen_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert process_body is not None, \
            "NORMALISE_CNV_QUALITY_SCORES process not found"
        assert "'.annotated.vcf.gz'" not in process_body and \
               '".annotated.vcf.gz"' not in process_body, (
            "NORMALISE_CNV_QUALITY_SCORES in DRAGEN module must NOT strip "
            ".annotated.vcf.gz; the input is now a .sorted.vcf.gz file"
        )


# ---------------------------------------------------------------------------
# 5. main.nf RUN_DRAGEN uses normalised_vcf for EVALUATE
# ---------------------------------------------------------------------------

class TestRunDragenEvaluateInput:
    """RUN_DRAGEN sub-workflow must pass normalised_vcf to EVALUATE."""

    def test_uses_normalised_vcf_for_evaluate(self, main_text):
        run_dragen = re.search(
            r'workflow RUN_DRAGEN\s*\{(.+?)(?=\nworkflow |\Z)',
            main_text, re.DOTALL
        )
        assert run_dragen is not None, "RUN_DRAGEN workflow not found in main.nf"
        body = run_dragen.group(1)
        assert 'DRAGEN.out.normalised_vcf' in body, (
            "RUN_DRAGEN must pass DRAGEN.out.normalised_vcf to EVALUATE "
            "so downstream evaluation uses QUAL_norm in the QUAL column"
        )

    def test_does_not_use_sorted_vcf_for_evaluate(self, main_text):
        run_dragen = re.search(
            r'workflow RUN_DRAGEN\s*\{(.+?)(?=\nworkflow |\Z)',
            main_text, re.DOTALL
        )
        assert run_dragen is not None, "RUN_DRAGEN workflow not found in main.nf"
        body = run_dragen.group(1)
        assert 'DRAGEN.out.sorted_vcf' not in body, (
            "RUN_DRAGEN must NOT use DRAGEN.out.sorted_vcf for EVALUATE; "
            "use DRAGEN.out.normalised_vcf to ensure QUAL_norm is used"
        )

    def test_does_not_use_annotated_vcfs_for_evaluate(self, main_text):
        run_dragen = re.search(
            r'workflow RUN_DRAGEN\s*\{(.+?)(?=\nworkflow |\Z)',
            main_text, re.DOTALL
        )
        assert run_dragen is not None, "RUN_DRAGEN workflow not found in main.nf"
        body = run_dragen.group(1)
        assert 'DRAGEN.out.annotated_vcfs' not in body, (
            "RUN_DRAGEN must NOT use DRAGEN.out.annotated_vcfs for EVALUATE; "
            "use DRAGEN.out.normalised_vcf instead"
        )


# ---------------------------------------------------------------------------
# 6. ADD_DRAGEN_TOOL_ANNOTATION strips .cnv.vcf.gz (DRAGEN Germline Enrichment)
# ---------------------------------------------------------------------------

class TestDragenAnnotationCnvVcfGzNaming:
    """ADD_DRAGEN_TOOL_ANNOTATION must handle ${sample_id}.cnv.vcf.gz inputs.

    DRAGEN Germline Enrichment produces CNV VCF files named
    ${sample_id}.cnv.vcf.gz.  The sample_name derivation in the process must
    strip '.cnv.vcf.gz' (and '.cnv.vcf') before constructing the output
    filename so that the annotated VCF is named correctly.
    """

    def test_strips_cnv_vcf_gz_in_sample_name(self, dragen_text):
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found in modules-icav2-dragen.nf"
        # In Nextflow script blocks backslashes are doubled; match either form.
        has_cnv = (
            r'\\.cnv\\.vcf\\.gz' in process_body   # escaped form inside script block
            or '.cnv.vcf.gz' in process_body        # unescaped form
        )
        assert has_cnv, (
            "ADD_DRAGEN_TOOL_ANNOTATION must strip .cnv.vcf.gz when deriving "
            "sample_name to support DRAGEN Germline Enrichment output naming"
        )

    def test_find_targets_only_cnv_vcf_gz(self, dragen_text):
        """find must target only *.cnv.vcf.gz, not all VCF files.

        The DRAGEN Germline Enrichment output directory contains many VCF files
        (clean_CNVs.vcf, filtered_output.vcf, filtered_output_20.vcf, etc.).
        Only the primary CNV VCF – ${sample_id}.cnv.vcf.gz – should be
        INFO-tagged and fed into the normalisation chain.
        """
        process_body = _extract_process(dragen_text, 'ADD_DRAGEN_TOOL_ANNOTATION')
        assert process_body is not None, \
            "ADD_DRAGEN_TOOL_ANNOTATION process not found in modules-icav2-dragen.nf"
        assert '*.cnv.vcf.gz' in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must use find -name '*.cnv.vcf.gz' to target "
            "only the primary CNV VCF produced by DRAGEN Germline Enrichment"
        )
        # Ensure the broad patterns that would pick up clean_CNVs.vcf etc. are gone
        assert '-name "*.vcf"' not in process_body and "-name '*.vcf'" not in process_body, (
            "ADD_DRAGEN_TOOL_ANNOTATION must NOT use find -name '*.vcf'; "
            "this would match clean_CNVs.vcf, filtered_output.vcf, etc."
        )


# ---------------------------------------------------------------------------
# 7. Full naming-convention chain: _DRAGEN.normalised.vcf.gz
# ---------------------------------------------------------------------------

class TestDragenNormalisedNamingConvention:
    """DRAGEN normalised VCF must follow the ${sample_id}_DRAGEN.normalised.vcf.gz convention.

    main.nf build_feature_inputs references normalised DRAGEN VCFs as
    '${sample_id}_DRAGEN.normalised.vcf.gz'.  The three-step chain in the
    DRAGEN module must produce exactly this name:

      ADD_DRAGEN_TOOL_ANNOTATION  →  ${sample_id}_DRAGEN.annotated.vcf.gz
      BGZIP_SORT_INDEX_VCF        →  ${sample_id}_DRAGEN.sorted.vcf.gz
      NORMALISE_CNV_QUALITY_SCORES→  ${sample_id}_DRAGEN.normalised.vcf.gz
    """

    def test_bgzip_strips_annotated_suffix_to_derive_sample_name(self, dragen_text):
        """BGZIP_SORT_INDEX_VCF must strip .annotated.vcf.gz (not .sorted.vcf.gz)."""
        process_body = _extract_process(dragen_text, 'BGZIP_SORT_INDEX_VCF')
        assert process_body is not None, "BGZIP_SORT_INDEX_VCF not found"
        assert "'.annotated.vcf.gz'" in process_body or '".annotated.vcf.gz"' in process_body, (
            "BGZIP_SORT_INDEX_VCF must strip .annotated.vcf.gz to derive sample_name, "
            "so that the sorted VCF is named ${sample_id}_DRAGEN.sorted.vcf.gz"
        )

    def test_main_nf_expects_dragen_normalised_vcf(self, main_text):
        """main.nf build_feature_inputs must reference ${sample_id}_DRAGEN.normalised.vcf.gz."""
        assert '_DRAGEN.normalised.vcf.gz' in main_text, (
            "main.nf build_feature_inputs must reference the DRAGEN normalised VCF as "
            "'${sample_id}_DRAGEN.normalised.vcf.gz'"
        )
        assert '_DRAGEN_output.normalised.vcf.gz' not in main_text, (
            "main.nf must NOT reference '_DRAGEN_output.normalised.vcf.gz'; "
            "the convention is '${sample_id}_DRAGEN.normalised.vcf.gz'"
        )
