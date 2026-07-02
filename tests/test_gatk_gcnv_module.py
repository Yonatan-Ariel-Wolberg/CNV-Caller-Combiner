#!/usr/bin/env python3
"""
Tests for the GATK-gCNV Nextflow module (modules/callers/modules-gatk-gcnv.nf).

Validates that the module follows the GATK gCNV cohort-calling best practices
for exome sequencing data:
  1. Exome-appropriate default parameters are centrally declared in nextflow.config
     (bin_length=0, padding=250, scatter_count=5000, is_wgs=false).
  2. DetermineGermlineContigPloidy uses --run-mode COHORT.
  3. GermlineCNVCaller uses per-shard unique output prefixes derived from the
     shard directory name (not the static "shard" prefix that caused staging
     conflicts when shards were collected).
  4. PostprocessGermlineCNVCalls receives shards sorted by name and uses Groovy
     GString interpolation (${call_args} / ${model_args}), not escaped shell
     variable references ($call_args / $model_args).
  5. BGZIP_SORT_INDEX_VCF accepts a (sample_id, vcf) tuple and produces output
     files named <sample_id>_GCNV_genotyped_segments.sorted.vcf.gz so that
     gather_vcfs() in main.nf can extract the bare sample ID.
  6. The GATK_GCNV workflow sorts counts by filename and sorted-lists the
     sample index channel so that sample indices are deterministically assigned.
"""

import os
import re
import pytest


NF_MODULE = "modules/callers/modules-gatk-gcnv.nf"
NEXTFLOW_CONFIG = "nextflow.config"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _read_module(repo_root="."):
    path = f"{repo_root}/{NF_MODULE}"
    with open(path) as fh:
        return fh.read()


def _read_config(repo_root="."):
    path = f"{repo_root}/{NEXTFLOW_CONFIG}"
    with open(path) as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def module_text():
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    return _read_module(repo_root)


@pytest.fixture(scope="module")
def global_config_text():
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    content = _read_config(repo_root)
    return content.split("profiles {")[0]


class TestParamDeclarationStrategy:
    """GATK-gCNV module should consume globally-declared params."""

    def test_module_has_no_strict_params_block(self, module_text):
        assert "params {" not in module_text, (
            "modules-gatk-gcnv.nf should not define a strict module-level params block; "
            "pipeline-facing params are centralized globally."
        )


# ===========================================================================
# 1. Default parameter values suitable for exome sequencing
# ===========================================================================

class TestDefaultParams:
    """Global params defaults must reflect exome best-practice defaults."""

    def test_bin_length_is_zero(self, global_config_text):
        """bin_length must be 0 for exome (no binning; each target is its own interval)."""
        match = re.search(r"bin_length\s*=\s*(\d+)", global_config_text)
        assert match, "bin_length parameter not found in global params block"
        assert match.group(1) == "0", (
            f"bin_length should be 0 for exome data, got {match.group(1)}"
        )

    def test_padding_is_250(self, global_config_text):
        """padding must be 250 bp per GATK gCNV exome best practice."""
        match = re.search(r"padding\s*=\s*(\d+)", global_config_text)
        assert match, "padding parameter not found in global params block"
        assert match.group(1) == "250", (
            f"padding should be 250 for exome data, got {match.group(1)}"
        )

    def test_scatter_count_is_5000(self, global_config_text):
        """scatter_count (SCATTER_CONTENT) must match the JSON params default of 5000."""
        match = re.search(r"scatter_count\s*=\s*(\d+)", global_config_text)
        assert match, "scatter_count parameter not found in global params block"
        assert match.group(1) == "5000", (
            f"scatter_count should be 5000 (intervals per shard), got {match.group(1)}"
        )

    def test_is_wgs_is_false(self, global_config_text):
        """is_wgs must be false (exome, not WGS)."""
        match = re.search(r"is_wgs\s*=\s*(\S+)", global_config_text)
        assert match, "is_wgs parameter not found in global params block"
        assert match.group(1) == "false", (
            f"is_wgs should be false for exome data, got {match.group(1)}"
        )


# ===========================================================================
# 2. DetermineGermlineContigPloidy must use --run-mode COHORT
# ===========================================================================

class TestDetermineGermlineContigPloidy:
    """DetermineGermlineContigPloidy must explicitly specify COHORT run mode."""

    def test_run_mode_cohort_present(self, module_text):
        """--run-mode COHORT must appear in the DetermineGermlineContigPloidy command."""
        # Find the DETERMINE_PLOIDY_COHORT process block
        proc_match = re.search(
            r"process DETERMINE_PLOIDY_COHORT\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "DETERMINE_PLOIDY_COHORT process not found"
        proc_body = proc_match.group(1)
        assert "--run-mode COHORT" in proc_body, (
            "DetermineGermlineContigPloidy must include --run-mode COHORT "
            "for cohort-mode CNV calling"
        )


# ===========================================================================
# 3. GermlineCNVCaller: unique per-shard output prefix
# ===========================================================================

class TestGermlineCNVCaller:
    """GermlineCNVCaller must produce uniquely-named shard directories."""

    def test_output_prefix_uses_shard_parent_name(self, module_text):
        """The --output-prefix must be derived from interval_shard.parent.name,
        not the static string 'shard' which caused all shards to collide."""
        proc_match = re.search(
            r"process GERMLINE_CNV_CALLER_COHORT\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "GERMLINE_CNV_CALLER_COHORT process not found"
        proc_body = proc_match.group(1)
        # The script must NOT use the static prefix "shard"
        assert "--output-prefix shard" not in proc_body, (
            "Static --output-prefix 'shard' causes all shard directories to have "
            "the same name. Use interval_shard.parent.name instead."
        )
        # The output declarations must NOT use the static names 'shard-calls'/'shard-model'
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", proc_body, re.DOTALL
        )
        assert output_section, "output section not found in GERMLINE_CNV_CALLER_COHORT"
        out_text = output_section.group(1)
        assert '"shard-calls"' not in out_text, (
            "Static output path 'shard-calls' causes naming conflicts; "
            "use a shard-specific prefix."
        )
        assert '"shard-model"' not in out_text, (
            "Static output path 'shard-model' causes naming conflicts; "
            "use a shard-specific prefix."
        )

    def test_output_uses_interval_shard_parent_name(self, module_text):
        """Output paths must reference interval_shard.parent.name for unique naming."""
        proc_match = re.search(
            r"process GERMLINE_CNV_CALLER_COHORT\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "GERMLINE_CNV_CALLER_COHORT process not found"
        proc_body = proc_match.group(1)
        assert "interval_shard.parent.name" in proc_body, (
            "Output prefix must be derived from interval_shard.parent.name "
            "to guarantee unique shard directory names."
        )


# ===========================================================================
# 4. PostprocessGermlineCNVCalls: correct Groovy interpolation and shard sorting
# ===========================================================================

class TestPostprocessCalls:
    """POSTPROCESS_CALLS must pass shard args via Groovy interpolation, not shell vars."""

    def test_no_escaped_dollar_call_args(self, module_text):
        r"""'\$call_args' in the script block is a shell variable (always empty).
        The correct form is '${call_args}' which uses Groovy GString interpolation."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert r"\$call_args" not in proc_body, (
            r"'\$call_args' is an unset shell variable. "
            "Use '${call_args}' for Groovy GString interpolation."
        )

    def test_no_escaped_dollar_model_args(self, module_text):
        r"""'\$model_args' in the script block is a shell variable (always empty)."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert r"\$model_args" not in proc_body, (
            r"'\$model_args' is an unset shell variable. "
            "Use '${model_args}' for Groovy GString interpolation."
        )

    def test_call_args_groovy_interpolation(self, module_text):
        """${call_args} must appear in the script heredoc for correct interpolation."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert "${call_args}" in proc_body, (
            "${call_args} (Groovy interpolation) not found in POSTPROCESS_CALLS script"
        )

    def test_model_args_groovy_interpolation(self, module_text):
        """${model_args} must appear in the script heredoc for correct interpolation."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert "${model_args}" in proc_body, (
            "${model_args} (Groovy interpolation) not found in POSTPROCESS_CALLS script"
        )

    def test_shards_are_sorted(self, module_text):
        """Shards must be sorted before building CLI args so genomic order is preserved."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert ".sort" in proc_body, (
            "Call/model shards must be sorted by name before building "
            "--calls-shard-path / --model-shard-path arguments so that "
            "PostprocessGermlineCNVCalls receives them in correct genomic order."
        )

    def test_collect_closure_no_escaped_it(self, module_text):
        r"""The collect closure must not use '\$it' which produces the literal
        string '$it' rather than the Groovy closure variable value."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        proc_body = proc_match.group(1)
        assert r"\$it" not in proc_body, (
            r"'\$it' in a collect closure produces the literal text '$it', not the "
            "path value. Use an explicit variable name like '{ c -> \"...${c}\" }'."
        )

    def test_output_emits_tuple_with_sample_id(self, module_text):
        """POSTPROCESS_CALLS final_vcf output must be a tuple(val(sample_id), path(...))
        so that BGZIP_SORT_INDEX_VCF can name the file using the _GCNV convention."""
        proc_match = re.search(
            r"process POSTPROCESS_CALLS\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "POSTPROCESS_CALLS process not found"
        # The output section
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)",
            proc_match.group(1),
            re.DOTALL,
        )
        assert output_section, "output section not found in POSTPROCESS_CALLS"
        out_text = output_section.group(1)
        # Must be a tuple emit for final_vcf
        assert re.search(r"tuple\s+val\(sample_id\)", out_text), (
            "POSTPROCESS_CALLS final_vcf must emit tuple(val(sample_id), path(...)) "
            "so BGZIP_SORT_INDEX_VCF receives the sample ID for correct output naming."
        )


# ===========================================================================
# 5. BGZIP_SORT_INDEX_VCF: correct input tuple and _GCNV naming
# ===========================================================================

class TestBgzipSortIndexVcf:
    """BGZIP_SORT_INDEX_VCF must accept a tuple and produce _GCNV-tagged output."""

    def test_input_is_tuple_with_sample_id(self, module_text):
        """Input must be tuple(val(sample_id), path(vcf_file)) to receive sample ID."""
        proc_match = re.search(
            r"process BGZIP_SORT_INDEX_VCF\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "BGZIP_SORT_INDEX_VCF process not found"
        proc_body = proc_match.group(1)
        input_section = re.search(
            r"input:\s*(.+?)(?=\s*output:)", proc_body, re.DOTALL
        )
        assert input_section, "input section not found in BGZIP_SORT_INDEX_VCF"
        in_text = input_section.group(1)
        assert re.search(r"tuple\s+val\(sample_id\)", in_text), (
            "BGZIP_SORT_INDEX_VCF input must be tuple(val(sample_id), path(vcf_file))"
        )

    def test_output_filename_contains_gcnv_tag(self, module_text):
        """Output filename must contain '_GCNV' so gather_vcfs() regex extracts
        the correct sample ID when building consensus VCF inputs."""
        proc_match = re.search(
            r"process BGZIP_SORT_INDEX_VCF\s*\{(.+?)(?=\nprocess |\nworkflow )",
            module_text,
            re.DOTALL,
        )
        assert proc_match, "BGZIP_SORT_INDEX_VCF process not found"
        proc_body = proc_match.group(1)
        assert "_GCNV" in proc_body, (
            "BGZIP_SORT_INDEX_VCF output filename must include '_GCNV' so that "
            "gather_vcfs() regex '_(CANOES|...|GCNV|...).*' strips it correctly "
            "and recovers the bare sample ID."
        )


# ===========================================================================
# 6. Workflow: deterministic sample ordering
# ===========================================================================

class TestWorkflowOrdering:
    """The GATK_GCNV workflow must sort counts and sample indices consistently."""

    def test_counts_channel_is_sorted(self, module_text):
        """Counts must be sorted by filename before being passed to
        FilterIntervals / DetermineGermlineContigPloidy / GermlineCNVCaller."""
        wf_match = re.search(
            r"workflow GATK_GCNV\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "GATK_GCNV workflow not found"
        wf_body = wf_match.group(1)
        # toSortedList with a comparator based on file name
        assert "toSortedList" in wf_body, (
            "The GATK_GCNV workflow must sort the counts channel with toSortedList "
            "so that FilterIntervals, DetermineGermlineContigPloidy, and "
            "GermlineCNVCaller all receive count files in the same deterministic order."
        )

    def test_indexed_samples_channel_is_sorted(self, module_text):
        """Sample indices must be assigned in alphabetical (sorted) order to match
        the order used when building counts_files.list for DetermineGermlineContigPloidy."""
        wf_match = re.search(
            r"workflow GATK_GCNV\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "GATK_GCNV workflow not found"
        wf_body = wf_match.group(1)
        # The indexed_samples_ch construction must use toSortedList (not just toList)
        assert "toSortedList" in wf_body, (
            "indexed_samples_ch must use toSortedList() so sample indices match "
            "the sorted order of count files passed to DetermineGermlineContigPloidy."
        )
        # Specifically for the indexed_samples_ch block, toList() must not be used
        idx_block = re.search(
            r"indexed_samples_ch\s*=\s*bam_ch(.+?)indexed_samples_ch\b",
            wf_body,
            re.DOTALL,
        )
        if idx_block:
            assert ".toList()" not in idx_block.group(1), (
                "Replace .toList() with .toSortedList() for the indexed samples channel "
                "to guarantee consistent sample indexing."
            )


# ===========================================================================
# 7. NORMALISE_CNV_QUALITY_SCORES: shared process from modules-common.nf
# ===========================================================================

class TestNormaliseProcess:
    """NORMALISE_CNV_QUALITY_SCORES must be imported from modules-common.nf,
    not redefined locally in the GATK module."""

    def test_normalise_process_not_defined_locally(self, module_text):
        """The GATK module must not contain a local definition of
        NORMALISE_CNV_QUALITY_SCORES; it should be included from modules-common.nf."""
        # Count how many times 'process NORMALISE_CNV_QUALITY_SCORES' appears
        local_defs = re.findall(r"^\s*process\s+NORMALISE_CNV_QUALITY_SCORES\b", module_text, re.MULTILINE)
        assert len(local_defs) == 0, (
            "NORMALISE_CNV_QUALITY_SCORES must not be defined locally in the GATK "
            "module; include it from modules-common.nf instead."
        )

    def test_normalise_included_from_common(self, module_text):
        """NORMALISE_CNV_QUALITY_SCORES must be included from modules-common.nf."""
        assert re.search(
            r"include\s*\{[^}]*NORMALISE_CNV_QUALITY_SCORES[^}]*\}.*modules-common\.nf",
            module_text,
        ), (
            "modules-gatk-gcnv.nf must include NORMALISE_CNV_QUALITY_SCORES "
            "from ../common/modules-common.nf"
        )

    def test_workflow_normalise_call_passes_caller_gatk(self, module_text):
        """The GATK_GCNV workflow must pass 'GATK' as the caller_name argument
        when invoking NORMALISE_CNV_QUALITY_SCORES (shared process signature)."""
        wf_match = re.search(
            r"workflow GATK_GCNV\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "GATK_GCNV workflow not found"
        wf_body = wf_match.group(1)
        # Find the line(s) containing the NORMALISE_CNV_QUALITY_SCORES call
        normalise_call = re.search(
            r"NORMALISE_CNV_QUALITY_SCORES\s*\(.+?'GATK'",
            wf_body,
            re.DOTALL,
        )
        assert normalise_call, (
            "NORMALISE_CNV_QUALITY_SCORES call must pass 'GATK' as the caller_name "
            "argument (second positional arg) to match the shared process signature."
        )

    def test_workflow_normalise_call_passes_dir_suffix(self, module_text):
        """The GATK_GCNV workflow must pass 'out_GCNV' as the dir_suffix argument
        when invoking NORMALISE_CNV_QUALITY_SCORES (shared process signature)."""
        wf_match = re.search(
            r"workflow GATK_GCNV\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "GATK_GCNV workflow not found"
        wf_body = wf_match.group(1)
        # Find the line(s) containing the NORMALISE_CNV_QUALITY_SCORES call
        normalise_call = re.search(
            r"NORMALISE_CNV_QUALITY_SCORES\s*\(.+?'out_GCNV'",
            wf_body,
            re.DOTALL,
        )
        assert normalise_call, (
            "NORMALISE_CNV_QUALITY_SCORES call must pass 'out_GCNV' as the dir_suffix "
            "argument (third positional arg) to match the shared process signature."
        )


# ===========================================================================
# 8. gcnv workflow numeric param validation in main.nf
# ===========================================================================

class TestGcnvParamValidation:
    """validate_required_params in main.nf must validate numeric parameters for
    the gcnv workflow: bin_length >= 0, padding >= 0, scatter_count > 0."""

    @pytest.fixture(scope="class")
    def main_text(self):
        import os
        main_nf = os.path.join(os.path.dirname(__file__), "..", "main.nf")
        with open(main_nf) as fh:
            return fh.read()

    def test_gcnv_validation_block_present(self, main_text):
        """validate_required_params must include a gcnv-specific validation block."""
        assert "workflow_name == 'gcnv'" in main_text, (
            "validate_required_params must contain a block checking "
            "workflow_name == 'gcnv' to apply gcnv-specific param validation"
        )

    def test_bin_length_non_negative_check(self, main_text):
        """bin_length must be validated as non-negative (>= 0) for gcnv workflow."""
        assert "bin_length_val < 0" in main_text or \
               "--bin_length must be a non-negative integer" in main_text, (
            "validate_required_params must check that bin_length >= 0 for gcnv workflow"
        )

    def test_bin_length_error_message(self, main_text):
        """The bin_length validation error must identify the param and the constraint."""
        assert "--bin_length must be a non-negative integer for --workflow gcnv" in main_text, (
            "Validation error message for bin_length must say "
            "'--bin_length must be a non-negative integer for --workflow gcnv'"
        )

    def test_padding_non_negative_check(self, main_text):
        """padding must be validated as non-negative (>= 0) for gcnv workflow."""
        assert "padding_val < 0" in main_text or \
               "--padding must be a non-negative integer" in main_text, (
            "validate_required_params must check that padding >= 0 for gcnv workflow"
        )

    def test_padding_error_message(self, main_text):
        """The padding validation error must identify the param and the constraint."""
        assert "--padding must be a non-negative integer for --workflow gcnv" in main_text, (
            "Validation error message for padding must say "
            "'--padding must be a non-negative integer for --workflow gcnv'"
        )

    def test_scatter_count_positive_check(self, main_text):
        """scatter_count must be validated as strictly positive (> 0) for gcnv workflow."""
        assert "scatter_count_val <= 0" in main_text or \
               "--scatter_count must be a positive integer" in main_text, (
            "validate_required_params must check that scatter_count > 0 for gcnv workflow"
        )

    def test_scatter_count_error_message(self, main_text):
        """The scatter_count validation error must identify the param and constraint."""
        assert "--scatter_count must be a positive integer for --workflow gcnv" in main_text, (
            "Validation error message for scatter_count must say "
            "'--scatter_count must be a positive integer for --workflow gcnv'"
        )
