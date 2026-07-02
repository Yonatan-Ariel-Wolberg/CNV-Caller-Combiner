#!/usr/bin/env python3
"""
Tests for the XHMM Nextflow module (modules/callers/modules-xhmm.nf).

Validates that the module correctly implements the XHMM CNV-calling pipeline
for exome sequencing from a cohort, producing per-sample VCF files that:
  1. Only contain CNVs actually called for that sample (GT-based filtering).
  2. Have the 'chr' prefix in the CHROM column of output VCFs.
  3. Use output filenames containing '_XHMM_' so gather_vcfs() can extract
     sample IDs (see main.nf get_id regex).
  4. Handle references both with and without chromosome 'chr' prefix
     in the CALC_GC_XHMM step.
"""

import os
import re

import pytest

NF_MODULE = "modules/callers/modules-xhmm.nf"


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


def _extract_script(process_body):
    """Return the heredoc script block content from a process body."""
    match = re.search(r'"""\s*(.+?)\s*"""', process_body, re.DOTALL)
    return match.group(1) if match else ""


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def module_text():
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    return _read_module(repo_root)


@pytest.fixture(scope="module")
def split_vcf_script(module_text):
    body = _extract_process(module_text, "SPLIT_VCF")
    assert body is not None, "SPLIT_VCF process not found in modules-xhmm.nf"
    return _extract_script(body)


@pytest.fixture(scope="module")
def calc_gc_script(module_text):
    body = _extract_process(module_text, "CALC_GC_XHMM")
    assert body is not None, "CALC_GC_XHMM process not found in modules-xhmm.nf"
    return _extract_script(body)


@pytest.fixture(scope="module")
def filter_cnvs_script(module_text):
    body = _extract_process(module_text, "FILTER_XHMM_CNVS")
    assert body is not None, "FILTER_XHMM_CNVS process not found in modules-xhmm.nf"
    return _extract_script(body)


# ---------------------------------------------------------------------------
# SPLIT_VCF tests
# ---------------------------------------------------------------------------

class TestSplitVcf:
    """SPLIT_VCF must produce per-sample VCFs named *_XHMM_output.vcf,
    filter by GT (not by FORMAT field missing-value checks), and add
    the 'chr' prefix to the CHROM column."""

    def test_output_glob_contains_xhmm(self, split_vcf_script):
        """Output file glob must contain '_XHMM_output.vcf' so gather_vcfs()
        can strip '_XHMM_...' to recover the sample ID."""
        assert "_XHMM_output.vcf" in split_vcf_script, (
            "SPLIT_VCF output files must be named *_XHMM_output.vcf "
            "(required for gather_vcfs() regex in main.nf)"
        )

    def test_gt_based_filter_used(self, split_vcf_script):
        """SPLIT_VCF must filter by GT (keep only non-ref, non-missing calls)
        rather than checking whether FORMAT float fields equal '.'."""
        # Correct approach: include only variants where GT != 0 (ref) and GT != . (missing)
        assert 'GT!="."' in split_vcf_script, (
            "SPLIT_VCF must exclude missing genotype (GT=\".\") records "
            "to keep only actual CNV calls"
        )
        assert 'GT!="0"' in split_vcf_script, (
            "SPLIT_VCF must exclude reference genotype (GT=\"0\") records "
            "so only samples with actual CNV calls are retained"
        )

    def test_no_format_field_string_comparison(self, split_vcf_script):
        """SPLIT_VCF must NOT use 'FORMAT/EQ=\".\"' style filtering, which is
        fragile for float FORMAT fields; GT-based filtering is correct."""
        assert 'FORMAT/EQ="."' not in split_vcf_script, (
            "SPLIT_VCF should not filter using FORMAT/EQ=\".\" string comparison; "
            "use GT-based filtering instead"
        )
        assert 'FORMAT/SQ="."' not in split_vcf_script, (
            "SPLIT_VCF should not filter using FORMAT/SQ=\".\" string comparison"
        )
        assert 'FORMAT/NDQ="."' not in split_vcf_script, (
            "SPLIT_VCF should not filter using FORMAT/NDQ=\".\" string comparison"
        )

    def test_chr_prefix_added(self, split_vcf_script):
        """SPLIT_VCF must add 'chr' prefix to the CHROM column of output VCFs
        so that downstream tools and the final VCF have consistent chromosome names."""
        assert "chr" in split_vcf_script and ("!~ /^chr/" in split_vcf_script or
                                               "!~/^chr/" in split_vcf_script), (
            "SPLIT_VCF must add 'chr' prefix to CHROM field when it is absent "
            "(e.g., via awk with $1 !~ /^chr/)"
        )

    def test_splits_by_sample_using_bcftools_query(self, split_vcf_script):
        """SPLIT_VCF must use 'bcftools query -l' to get the list of samples
        from the multi-sample VCF."""
        assert "bcftools query -l" in split_vcf_script, (
            "SPLIT_VCF must use 'bcftools query -l' to enumerate samples "
            "from the multi-sample XHMM VCF"
        )

    def test_uses_bcftools_view_per_sample(self, split_vcf_script):
        """SPLIT_VCF must call 'bcftools view -s' to extract each sample's column."""
        assert "bcftools view -s" in split_vcf_script, (
            "SPLIT_VCF must use 'bcftools view -s <sample>' to subset the VCF "
            "to a single sample before filtering"
        )


# ---------------------------------------------------------------------------
# CALC_GC_XHMM tests
# ---------------------------------------------------------------------------

class TestCalcGcXhmm:
    """CALC_GC_XHMM must handle references with and without 'chr' prefix."""

    def test_does_not_use_grep_chr_only(self, calc_gc_script):
        """CALC_GC_XHMM must not rely solely on 'grep ^chr', which would
        silently produce empty output for non-chr-prefixed references."""
        # The old broken pattern was: grep "^chr" gc.tmp
        assert 'grep "^chr" gc.tmp' not in calc_gc_script, (
            "CALC_GC_XHMM must not use 'grep \"^chr\" gc.tmp' which fails "
            "for references without a 'chr' chromosome prefix; use "
            "'grep -v \"^@\"' or equivalent instead"
        )

    def test_skips_at_headers(self, calc_gc_script):
        """CALC_GC_XHMM must skip SAM-style '@'-prefixed header lines from
        GATK AnnotateIntervals output."""
        assert "grep -v" in calc_gc_script and '"^@"' in calc_gc_script, (
            "CALC_GC_XHMM must skip '@'-prefixed header lines in gc.tmp "
            "(e.g., grep -v '^@')"
        )

    def test_skips_contig_header(self, calc_gc_script):
        """CALC_GC_XHMM must skip the 'CONTIG' column-header line that GATK
        AnnotateIntervals writes at the top of the data block."""
        assert "CONTIG" in calc_gc_script, (
            "CALC_GC_XHMM must filter out the 'CONTIG' column-header line "
            "produced by GATK AnnotateIntervals"
        )

    def test_outputs_locus_gc_txt(self, calc_gc_script):
        """CALC_GC_XHMM must write DATA.locus_GC.txt."""
        assert "DATA.locus_GC.txt" in calc_gc_script, (
            "CALC_GC_XHMM must output DATA.locus_GC.txt"
        )

    def test_outputs_extreme_gc_targets(self, calc_gc_script):
        """CALC_GC_XHMM must write extreme_gc_targets.txt."""
        assert "extreme_gc_targets.txt" in calc_gc_script, (
            "CALC_GC_XHMM must output extreme_gc_targets.txt"
        )


# ---------------------------------------------------------------------------
# FILTER_XHMM_CNVS tests
# ---------------------------------------------------------------------------

class TestFilterXhmmCnvs:
    """FILTER_XHMM_CNVS must apply quality-score thresholds correctly."""

    def test_uses_bcftools_filter(self, filter_cnvs_script):
        """FILTER_XHMM_CNVS must call bcftools filter to apply quality thresholds."""
        assert "bcftools filter" in filter_cnvs_script, (
            "FILTER_XHMM_CNVS must use bcftools filter to apply quality thresholds"
        )

    def test_filters_on_eq_sq_ndq(self, filter_cnvs_script):
        """FILTER_XHMM_CNVS must reference EQ, SQ, and NDQ FORMAT fields."""
        assert "EQ" in filter_cnvs_script, "FILTER_XHMM_CNVS must check FORMAT/EQ"
        assert "SQ" in filter_cnvs_script, "FILTER_XHMM_CNVS must check FORMAT/SQ"
        assert "NDQ" in filter_cnvs_script, "FILTER_XHMM_CNVS must check FORMAT/NDQ"

    def test_output_naming_contains_xhmm(self, module_text):
        """Output glob of FILTER_XHMM_CNVS must produce filenames that still
        contain '_XHMM_' so that gather_vcfs() regex in main.nf works."""
        body = _extract_process(module_text, "FILTER_XHMM_CNVS")
        assert body is not None, "FILTER_XHMM_CNVS process not found"
        # The output glob is *_filtered.vcf; since input is *_XHMM_output.vcf,
        # output will be *_XHMM_output_filtered.vcf which contains _XHMM_
        script = _extract_script(body)
        # Check that the rename pattern preserves the base name
        assert "_filtered.vcf" in script, (
            "FILTER_XHMM_CNVS must name output files *_filtered.vcf "
            "so that '_XHMM_' remains in the filename from the input"
        )


# ---------------------------------------------------------------------------
# Workflow structure tests
# ---------------------------------------------------------------------------

class TestXhmmWorkflow:
    """The XHMM workflow block must wire processes in the correct order."""

    def _get_workflow_block(self, module_text):
        match = re.search(r'workflow XHMM \{(.+)', module_text, re.DOTALL)
        assert match, "workflow XHMM not found in modules-xhmm.nf"
        return match.group(1)

    def test_split_vcf_in_workflow(self, module_text):
        """SPLIT_VCF must be called in the XHMM workflow."""
        wf = self._get_workflow_block(module_text)
        assert "SPLIT_VCF" in wf, "SPLIT_VCF not called in workflow XHMM"

    def test_filter_xhmm_cnvs_after_split_vcf(self, module_text):
        """FILTER_XHMM_CNVS must appear after SPLIT_VCF in the workflow."""
        wf = self._get_workflow_block(module_text)
        split_pos = wf.find("SPLIT_VCF")
        filter_pos = wf.find("FILTER_XHMM_CNVS")
        assert split_pos != -1, "SPLIT_VCF not found in workflow XHMM"
        assert filter_pos != -1, "FILTER_XHMM_CNVS not found in workflow XHMM"
        assert split_pos < filter_pos, (
            "SPLIT_VCF must precede FILTER_XHMM_CNVS in the XHMM workflow"
        )

    def test_bgzip_sort_index_after_filter(self, module_text):
        """BGZIP_SORT_INDEX_VCF must appear after FILTER_XHMM_CNVS in workflow."""
        wf = self._get_workflow_block(module_text)
        filter_pos = wf.find("FILTER_XHMM_CNVS")
        bgzip_pos = wf.find("BGZIP_SORT_INDEX_VCF")
        assert filter_pos != -1, "FILTER_XHMM_CNVS not found in workflow XHMM"
        assert bgzip_pos != -1, "BGZIP_SORT_INDEX_VCF not found in workflow XHMM"
        assert filter_pos < bgzip_pos, (
            "FILTER_XHMM_CNVS must precede BGZIP_SORT_INDEX_VCF in XHMM workflow"
        )

    def test_emit_sorted_vcf(self, module_text):
        """XHMM workflow must emit sorted_vcf and sorted_vcf_index channels."""
        wf = self._get_workflow_block(module_text)
        assert "sorted_vcf" in wf, (
            "XHMM workflow emit block must include sorted_vcf"
        )
        assert "sorted_vcf_index" in wf, (
            "XHMM workflow emit block must include sorted_vcf_index"
        )
