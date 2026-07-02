#!/usr/bin/env python3
"""
Tests for the Truvari Nextflow module (modules/sv-mergers/modules-truvari.nf).

Validates that the module correctly implements a two-step VCF consensus
pipeline for exome CNV callers:

  1. MERGE_VCFS: sorts, bgzips, indexes, and concatenates per-caller VCFs
     into a single merged VCF (merged_<sample_id>.vcf.gz) with a .tbi index.

  2. COLLAPSE_VCFS: runs `truvari collapse` with the recommended parameters
     for exome CNV consensus (--intra, pctseq 0, pctsize 0.5, pctovl 0.5,
     sizemin 50, sizemax 2000000) and produces two output VCFs:
       - <sample_id>_truvari_merged.vcf   (remaining variants after collapsing)
       - <sample_id>_truvari_collapsed.vcf (collapsed/representative calls)

  3. TRUVARI workflow: chains MERGE_VCFS → COLLAPSE_VCFS and emits both
     'merged_vcf' and 'collapsed_vcf' output channels.

  4. COLLAPSE_VCFS uses publishDir to write outputs to out_TRUVARI/<sample_id>.

  5. Both processes use label 'truvari' to receive the correct container
     (docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0).
"""

import os
import re

import pytest


NF_MODULE = "modules/sv-mergers/modules-truvari.nf"


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


@pytest.fixture(scope="module")
def merge_vcfs_body(module_text):
    body = _extract_process(module_text, "MERGE_VCFS")
    assert body is not None, "MERGE_VCFS process not found in module"
    return body


@pytest.fixture(scope="module")
def collapse_vcfs_body(module_text):
    body = _extract_process(module_text, "COLLAPSE_VCFS")
    assert body is not None, "COLLAPSE_VCFS process not found in module"
    return body


@pytest.fixture(scope="module")
def truvari_workflow_body(module_text):
    match = re.search(
        r"workflow TRUVARI\s*\{(.+)",
        module_text,
        re.DOTALL,
    )
    assert match, "TRUVARI workflow not found in module"
    return match.group(1)


# ===========================================================================
# 1. MERGE_VCFS: sorting, bgzipping, indexing, and concatenating VCFs
# ===========================================================================

class TestMergeVcfsProcess:
    """MERGE_VCFS must sort, bgzip, index, and concatenate per-caller VCFs."""

    def test_process_exists(self, module_text):
        """MERGE_VCFS process must be defined in the module."""
        assert "process MERGE_VCFS" in module_text, (
            "MERGE_VCFS process not found in modules-truvari.nf"
        )

    def test_label_is_truvari(self, merge_vcfs_body):
        """MERGE_VCFS must use label 'truvari' to receive the correct container."""
        assert "label 'truvari'" in merge_vcfs_body, (
            "MERGE_VCFS must have label 'truvari' to use the Truvari container"
        )

    def test_input_is_sample_vcfs_tuple(self, merge_vcfs_body):
        """Input must be tuple(val(sample_id), path(vcfs)) for groupTuple output."""
        input_section = re.search(
            r"input:\s*(.+?)(?=\s*output:)", merge_vcfs_body, re.DOTALL
        )
        assert input_section, "input section not found in MERGE_VCFS"
        in_text = input_section.group(1)
        assert re.search(r"tuple\s+val\(sample_id\)", in_text), (
            "MERGE_VCFS input must start with val(sample_id)"
        )
        assert "path(vcfs)" in in_text, (
            "MERGE_VCFS input must include path(vcfs) for the per-caller VCF list"
        )

    def test_output_merged_vcf_gz(self, merge_vcfs_body):
        """Output must include a bgzipped merged VCF named merged_<sample_id>.vcf.gz."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", merge_vcfs_body, re.DOTALL
        )
        assert output_section, "output section not found in MERGE_VCFS"
        out_text = output_section.group(1)
        assert "merged_" in out_text and ".vcf.gz" in out_text, (
            "MERGE_VCFS output must include merged_<sample_id>.vcf.gz"
        )

    def test_output_includes_tbi_index(self, merge_vcfs_body):
        """Output must include a .tbi index for the merged VCF."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", merge_vcfs_body, re.DOTALL
        )
        assert output_section, "output section not found in MERGE_VCFS"
        out_text = output_section.group(1)
        assert ".vcf.gz.tbi" in out_text, (
            "MERGE_VCFS output must include the .tbi index file "
            "so COLLAPSE_VCFS can open the indexed merged VCF"
        )

    def test_output_emits_merged_data(self, merge_vcfs_body):
        """Output tuple must use emit: merged_data for downstream consumption."""
        assert "emit: merged_data" in merge_vcfs_body, (
            "MERGE_VCFS output must use 'emit: merged_data' "
            "so COLLAPSE_VCFS can reference MERGE_VCFS.out.merged_data"
        )

    def test_bgzip_used_to_compress(self, merge_vcfs_body):
        """Script must use bgzip to compress VCFs (required for tabix indexing)."""
        script_section = re.search(r'"""(.+?)"""', merge_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in MERGE_VCFS"
        assert "bgzip" in script_section.group(1), (
            "MERGE_VCFS script must use bgzip to compress VCFs "
            "before sorting and indexing with tabix"
        )

    def test_bcftools_sort_used(self, merge_vcfs_body):
        """Script must sort each VCF with bcftools sort before concatenation."""
        script_section = re.search(r'"""(.+?)"""', merge_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in MERGE_VCFS"
        assert "bcftools sort" in script_section.group(1), (
            "MERGE_VCFS script must use bcftools sort to sort each per-caller VCF "
            "so that tabix can index it and bcftools concat receives sorted input"
        )

    def test_tabix_indexing_used(self, merge_vcfs_body):
        """Script must create .tbi index with tabix -p vcf."""
        script_section = re.search(r'"""(.+?)"""', merge_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in MERGE_VCFS"
        script_text = script_section.group(1)
        assert "tabix" in script_text, (
            "MERGE_VCFS script must use tabix to index the sorted VCFs"
        )
        assert "-p vcf" in script_text, (
            "tabix must be called with '-p vcf' preset for VCF files"
        )

    def test_script_uses_bcftools_concat(self, merge_vcfs_body):
        """Script must concatenate all sorted per-caller VCFs with bcftools concat."""
        script_section = re.search(r'"""(.+?)"""', merge_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in MERGE_VCFS"
        assert "bcftools concat" in script_section.group(1), (
            "MERGE_VCFS script must use bcftools concat to merge all per-caller VCFs "
            "into a single file before Truvari collapse"
        )

    def test_bcftools_concat_allows_overlaps(self, merge_vcfs_body):
        """bcftools concat must use -a flag to allow overlapping records from different callers."""
        script_section = re.search(r'"""(.+?)"""', merge_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in MERGE_VCFS"
        assert "bcftools concat -a" in script_section.group(1), (
            "bcftools concat must use the -a (--allow-overlaps) flag so records "
            "from different callers covering the same locus are retained"
        )


# ===========================================================================
# 2. COLLAPSE_VCFS: truvari collapse parameters
# ===========================================================================

class TestCollapseVcfsProcess:
    """COLLAPSE_VCFS must call truvari collapse with the correct parameters."""

    def test_process_exists(self, module_text):
        """COLLAPSE_VCFS process must be defined in the module."""
        assert "process COLLAPSE_VCFS" in module_text, (
            "COLLAPSE_VCFS process not found in modules-truvari.nf"
        )

    def test_label_is_truvari(self, collapse_vcfs_body):
        """COLLAPSE_VCFS must use label 'truvari' to receive the correct container."""
        assert "label 'truvari'" in collapse_vcfs_body, (
            "COLLAPSE_VCFS must have label 'truvari' to use the Truvari container"
        )

    def test_publishdir_points_to_out_truvari(self, collapse_vcfs_body):
        """publishDir must write to out_TRUVARI/<sample_id> for result organisation."""
        assert "out_TRUVARI" in collapse_vcfs_body, (
            "COLLAPSE_VCFS publishDir must point to out_TRUVARI/<sample_id> "
            "so results are organised alongside other caller outputs"
        )

    def test_input_accepts_merged_vcf_and_tbi(self, collapse_vcfs_body):
        """Input must accept (sample_id, merged_vcf, merged_vcf_tbi) from MERGE_VCFS."""
        input_section = re.search(
            r"input:\s*(.+?)(?=\s*output:)", collapse_vcfs_body, re.DOTALL
        )
        assert input_section, "input section not found in COLLAPSE_VCFS"
        in_text = input_section.group(1)
        assert re.search(r"tuple\s+val\(sample_id\)", in_text), (
            "COLLAPSE_VCFS input must start with val(sample_id)"
        )
        assert "merged_vcf" in in_text, (
            "COLLAPSE_VCFS input must include the merged VCF path"
        )

    def test_output_truvari_merged_vcf(self, collapse_vcfs_body):
        """Output must include <sample_id>_truvari_merged.vcf (remaining variants)."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", collapse_vcfs_body, re.DOTALL
        )
        assert output_section, "output section not found in COLLAPSE_VCFS"
        out_text = output_section.group(1)
        assert "_truvari_merged.vcf" in out_text, (
            "COLLAPSE_VCFS output must include <sample_id>_truvari_merged.vcf "
            "containing variants not collapsed by Truvari"
        )

    def test_output_truvari_collapsed_vcf(self, collapse_vcfs_body):
        """Output must include <sample_id>_truvari_collapsed.vcf (representative calls)."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", collapse_vcfs_body, re.DOTALL
        )
        assert output_section, "output section not found in COLLAPSE_VCFS"
        out_text = output_section.group(1)
        assert "_truvari_collapsed.vcf" in out_text, (
            "COLLAPSE_VCFS output must include <sample_id>_truvari_collapsed.vcf "
            "containing the representative (collapsed) calls from Truvari"
        )

    def test_output_emits_merged_vcf(self, collapse_vcfs_body):
        """Output must use emit: merged_vcf for the un-collapsed variants."""
        assert "emit: merged_vcf" in collapse_vcfs_body, (
            "COLLAPSE_VCFS must emit 'merged_vcf' so the TRUVARI workflow "
            "can expose it as an output channel"
        )

    def test_output_emits_collapsed_vcf(self, collapse_vcfs_body):
        """Output must use emit: collapsed_vcf for the representative calls."""
        assert "emit: collapsed_vcf" in collapse_vcfs_body, (
            "COLLAPSE_VCFS must emit 'collapsed_vcf' so the TRUVARI workflow "
            "can expose it as an output channel"
        )

    def test_truvari_collapse_command_present(self, collapse_vcfs_body):
        """Script must call 'truvari collapse'."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "truvari collapse" in script_section.group(1), (
            "COLLAPSE_VCFS script must call 'truvari collapse' to produce consensus calls"
        )

    def test_truvari_intra_flag(self, collapse_vcfs_body):
        """truvari collapse must use --intra to collapse within the same sample."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "--intra" in script_section.group(1), (
            "truvari collapse must include --intra to collapse variants "
            "from different callers that overlap within the same sample"
        )

    def test_truvari_pctseq_is_zero(self, collapse_vcfs_body):
        """truvari collapse --pctseq must be 0 (disable sequence similarity for CNVs)."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        script_text = script_section.group(1)
        match = re.search(r"--pctseq\s+(\S+)", script_text)
        assert match, "--pctseq not found in truvari collapse command"
        assert match.group(1) == "0", (
            f"--pctseq must be 0 for CNV calls (no sequence similarity required), "
            f"got {match.group(1)}"
        )

    def test_truvari_pctsize_is_half(self, collapse_vcfs_body):
        """truvari collapse --pctsize must be 0.5 (50% size overlap required)."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        script_text = script_section.group(1)
        match = re.search(r"--pctsize\s+(\S+)", script_text)
        assert match, "--pctsize not found in truvari collapse command"
        assert match.group(1) == "0.5", (
            f"--pctsize must be 0.5 (50% size similarity required to merge CNVs), "
            f"got {match.group(1)}"
        )

    def test_truvari_pctovl_is_half(self, collapse_vcfs_body):
        """truvari collapse --pctovl must be 0.5 (50% reciprocal overlap required)."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        script_text = script_section.group(1)
        match = re.search(r"--pctovl\s+(\S+)", script_text)
        assert match, "--pctovl not found in truvari collapse command"
        assert match.group(1) == "0.5", (
            f"--pctovl must be 0.5 (50% reciprocal overlap to merge CNVs), "
            f"got {match.group(1)}"
        )

    def test_truvari_sizemin_is_50(self, collapse_vcfs_body):
        """truvari collapse --sizemin must be 50 bp (ignore tiny SVs below exon resolution)."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        script_text = script_section.group(1)
        match = re.search(r"--sizemin\s+(\d+)", script_text)
        assert match, "--sizemin not found in truvari collapse command"
        assert int(match.group(1)) == 50, (
            f"--sizemin must be 50 bp, got {match.group(1)}"
        )

    def test_truvari_sizemax_is_2000000(self, collapse_vcfs_body):
        """truvari collapse --sizemax must be 2000000 bp (capture whole-chromosome CNVs)."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        script_text = script_section.group(1)
        match = re.search(r"--sizemax\s+(\d+)", script_text)
        assert match, "--sizemax not found in truvari collapse command"
        assert int(match.group(1)) == 2000000, (
            f"--sizemax must be 2000000 bp (2 Mb) to capture large chromosomal CNVs, "
            f"got {match.group(1)}"
        )

    def test_truvari_fast_cluster_flag(self, collapse_vcfs_body):
        """truvari collapse must use --fast-cluster for efficient large-cohort processing."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "--fast-cluster" in script_section.group(1), (
            "truvari collapse must include --fast-cluster for efficient processing "
            "of large multi-caller VCF files"
        )

    def test_truvari_input_flag(self, collapse_vcfs_body):
        """truvari collapse must use -i to specify the input merged VCF."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "-i " in script_section.group(1), (
            "truvari collapse must use -i to specify the merged input VCF"
        )

    def test_truvari_output_flag(self, collapse_vcfs_body):
        """truvari collapse must use -o to specify the merged output VCF."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "-o " in script_section.group(1), (
            "truvari collapse must use -o to specify the merged/remaining output VCF"
        )

    def test_truvari_collapsed_flag(self, collapse_vcfs_body):
        """truvari collapse must use -c to specify the collapsed representative VCF."""
        script_section = re.search(r'"""(.+?)"""', collapse_vcfs_body, re.DOTALL)
        assert script_section, "script block not found in COLLAPSE_VCFS"
        assert "-c " in script_section.group(1), (
            "truvari collapse must use -c to specify the collapsed representative VCF"
        )


# ===========================================================================
# 3. TRUVARI workflow: correct chaining and emit names
# ===========================================================================

class TestTruvariWorkflow:
    """TRUVARI workflow must chain MERGE_VCFS → COLLAPSE_VCFS and emit both outputs."""

    def test_workflow_exists(self, module_text):
        """TRUVARI workflow must be defined in the module."""
        assert "workflow TRUVARI" in module_text, (
            "TRUVARI workflow not found in modules-truvari.nf"
        )

    def test_workflow_takes_grouped_vcfs(self, truvari_workflow_body):
        """TRUVARI workflow take block must accept grouped_vcfs channel."""
        assert "grouped_vcfs" in truvari_workflow_body, (
            "TRUVARI workflow must declare 'grouped_vcfs' in its take block "
            "to receive (sample_id, [vcf1, vcf2, ...]) tuples from main.nf"
        )

    def test_merge_vcfs_called_in_workflow(self, truvari_workflow_body):
        """TRUVARI workflow must call MERGE_VCFS as first step."""
        assert "MERGE_VCFS" in truvari_workflow_body, (
            "TRUVARI workflow must call MERGE_VCFS to concatenate per-caller VCFs "
            "before collapsing"
        )

    def test_collapse_vcfs_called_in_workflow(self, truvari_workflow_body):
        """TRUVARI workflow must call COLLAPSE_VCFS after MERGE_VCFS."""
        assert "COLLAPSE_VCFS" in truvari_workflow_body, (
            "TRUVARI workflow must call COLLAPSE_VCFS to run truvari collapse "
            "on the merged VCF"
        )

    def test_collapse_receives_merge_output(self, truvari_workflow_body):
        """COLLAPSE_VCFS must receive MERGE_VCFS.out.merged_data as input."""
        assert "MERGE_VCFS.out.merged_data" in truvari_workflow_body, (
            "COLLAPSE_VCFS must be called with MERGE_VCFS.out.merged_data "
            "to chain the two processes correctly"
        )

    def test_merge_precedes_collapse_in_workflow(self, truvari_workflow_body):
        """MERGE_VCFS must appear before COLLAPSE_VCFS in the workflow body."""
        merge_pos = truvari_workflow_body.find("MERGE_VCFS")
        collapse_pos = truvari_workflow_body.find("COLLAPSE_VCFS")
        assert merge_pos != -1 and collapse_pos != -1, (
            "Both MERGE_VCFS and COLLAPSE_VCFS must appear in the TRUVARI workflow"
        )
        assert merge_pos < collapse_pos, (
            "MERGE_VCFS must be called before COLLAPSE_VCFS in the TRUVARI workflow"
        )

    def test_emit_merged_vcf(self, truvari_workflow_body):
        """TRUVARI workflow must emit 'merged_vcf' channel."""
        assert "merged_vcf" in truvari_workflow_body, (
            "TRUVARI workflow must emit 'merged_vcf' "
            "(the remaining variants after Truvari collapse)"
        )

    def test_emit_collapsed_vcf(self, truvari_workflow_body):
        """TRUVARI workflow must emit 'collapsed_vcf' channel."""
        assert "collapsed_vcf" in truvari_workflow_body, (
            "TRUVARI workflow must emit 'collapsed_vcf' "
            "(the representative calls from Truvari collapse)"
        )

    def test_merged_vcf_from_collapse_output(self, truvari_workflow_body):
        """merged_vcf emit must reference COLLAPSE_VCFS.out.merged_vcf."""
        assert "COLLAPSE_VCFS.out.merged_vcf" in truvari_workflow_body, (
            "TRUVARI workflow merged_vcf emit must be assigned from "
            "COLLAPSE_VCFS.out.merged_vcf"
        )

    def test_collapsed_vcf_from_collapse_output(self, truvari_workflow_body):
        """collapsed_vcf emit must reference COLLAPSE_VCFS.out.collapsed_vcf."""
        assert "COLLAPSE_VCFS.out.collapsed_vcf" in truvari_workflow_body, (
            "TRUVARI workflow collapsed_vcf emit must be assigned from "
            "COLLAPSE_VCFS.out.collapsed_vcf"
        )
