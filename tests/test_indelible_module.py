#!/usr/bin/env python3
"""
Tests for the INDELIBLE Nextflow module (modules/callers/modules-indelible.nf).

Validates the module against the reference implementations:
  * HurlesGroupSanger/indelible  (indelible.py subcommand CLI)
  * phelelani/nf-exomecnv        (Nextflow process structure)

Key correctness properties checked:

  1.  RUN_FETCH uses the correct indelible.py fetch flags and exports
      REF_PATH / REF_CACHE so that CRAM decoding works at runtime.

  2.  RUN_AGGREGATE passes all four required flags (--i --b --o --r).

  3.  RUN_SCORE passes the three required flags (--i --o --config).

  4.  RUN_DATABASE collects *.scored files via 'ls *.scored > scores.txt'
      and passes --priors and --tb to support cohort-level AF databases.

  5.  RUN_ANNOTATE uses the --d flag (database) and the correct I/O flags.

  6.  The INDELIBLE workflow broadcasts the single database channel to
      every sample by calling RUN_ANNOTATE with .first() on the database
      output.  Without .first(), Nextflow's default zip behaviour would
      only annotate the first sample.

  7.  RUN_DENOVO_TRIO uses --c (child annotation), --m (maternal BAM),
      and --p (paternal BAM).

  8.  RUN_DENOVO_MOM uses --c and --m only (no --p).

  9.  RUN_DENOVO_DAD uses --c and --p only (no --m).

  10. FILTER_INDELIBLE applies the canonical awk filter ($39 < 2 && $40 < 2)
      to the annotated TSV and preserves the sample ID in its output tuple.

  11. BGZIP_SORT_INDEX_VCF does NOT run a second bcftools-annotate pass to
      add TOOL=INDELIBLE — the converter already writes that field, so a
      second annotation would insert a duplicate ##INFO header line.

  12. BGZIP_SORT_INDEX_VCF uses bcftools sort with bgzip output (-O z) and
      tabix for indexing, producing *.sorted.vcf.gz and *.sorted.vcf.gz.tbi.

  13. The INDELIBLE workflow emits: filtered_cnvs, indelible_denovo,
      indelible_denovo_mom, indelible_denovo_dad, sorted_vcf, sorted_vcf_index.
"""

import re
import os
import pytest

NF_MODULE = "modules/callers/modules-indelible.nf"

# Repository root resolved once at import time
_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")


# ---------------------------------------------------------------------------
# Helper / fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def module_text():
    path = os.path.join(_REPO_ROOT, NF_MODULE)
    with open(path) as fh:
        return fh.read()


def _process_body(module_text, process_name):
    """Return the body of the named process block (everything after the opening brace).

    The lookahead ``(?=\\nprocess |\\nworkflow |\\Z)`` terminates the match at
    the next top-level ``process`` or ``workflow`` keyword, or at end-of-file.
    This avoids consuming text that belongs to the following declaration.
    """
    pattern = rf"process {re.escape(process_name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)"
    m = re.search(pattern, module_text, re.DOTALL)
    assert m, f"Process '{process_name}' not found in {NF_MODULE}"
    return m.group(1)


def _workflow_body(module_text, workflow_name="INDELIBLE"):
    """Return the body of the named workflow block."""
    pattern = rf"workflow {re.escape(workflow_name)}\s*\{{(.+)"
    m = re.search(pattern, module_text, re.DOTALL)
    assert m, f"Workflow '{workflow_name}' not found in {NF_MODULE}"
    return m.group(1)


# ===========================================================================
# 1. RUN_FETCH
# ===========================================================================

class TestRunFetch:
    """indelible.py fetch must receive --config, --i, --o and export CRAM env vars."""

    def test_ref_path_export(self, module_text):
        """REF_PATH must be exported so CRAM files can be decoded."""
        body = _process_body(module_text, "RUN_FETCH")
        assert "export REF_PATH=" in body, (
            "RUN_FETCH must export REF_PATH for CRAM reference decoding"
        )

    def test_ref_cache_export(self, module_text):
        """REF_CACHE must be exported alongside REF_PATH."""
        body = _process_body(module_text, "RUN_FETCH")
        assert "export REF_CACHE=" in body, (
            "RUN_FETCH must export REF_CACHE for CRAM reference decoding"
        )

    def test_fetch_command_flags(self, module_text):
        """indelible.py fetch must include --config, --i, and --o flags."""
        body = _process_body(module_text, "RUN_FETCH")
        assert "indelible.py fetch" in body
        assert "--config" in body
        assert " --i " in body
        assert " --o " in body

    def test_sc_reads_output(self, module_text):
        """RUN_FETCH must emit sc_reads so RUN_AGGREGATE can consume them."""
        body = _process_body(module_text, "RUN_FETCH")
        assert "emit: sc_reads" in body


# ===========================================================================
# 2. RUN_AGGREGATE
# ===========================================================================

class TestRunAggregate:
    """indelible.py aggregate requires --i, --b, --o, --r."""

    def test_aggregate_command_flags(self, module_text):
        body = _process_body(module_text, "RUN_AGGREGATE")
        assert "indelible.py aggregate" in body
        assert " --i " in body
        assert " --b " in body
        assert " --o " in body
        assert " --r " in body

    def test_counts_output(self, module_text):
        body = _process_body(module_text, "RUN_AGGREGATE")
        assert "emit: counts" in body


# ===========================================================================
# 3. RUN_SCORE
# ===========================================================================

class TestRunScore:
    """indelible.py score requires --i, --o, --config."""

    def test_score_command_flags(self, module_text):
        body = _process_body(module_text, "RUN_SCORE")
        assert "indelible.py score" in body
        assert " --i " in body
        assert " --o " in body
        assert "--config" in body

    def test_scores_and_database_in_outputs(self, module_text):
        """RUN_SCORE must emit both 'scores' (for annotate) and 'database_in' (for collect)."""
        body = _process_body(module_text, "RUN_SCORE")
        assert "emit: scores" in body
        assert "emit: database_in" in body


# ===========================================================================
# 4. RUN_DATABASE
# ===========================================================================

class TestRunDatabase:
    """indelible.py database requires --f, --o, --r, --priors, --config, --tb."""

    def test_database_command_flags(self, module_text):
        body = _process_body(module_text, "RUN_DATABASE")
        assert "indelible.py database" in body
        assert " --f " in body
        assert " --o " in body
        assert " --r " in body
        assert "--priors" in body
        assert "--config" in body
        assert "--tb" in body

    def test_ls_scored_files(self, module_text):
        """The database step must list *.scored files into a file-of-filenames."""
        body = _process_body(module_text, "RUN_DATABASE")
        assert "ls *.scored" in body, (
            "RUN_DATABASE must collect scored files with 'ls *.scored > scores.txt' "
            "so indelible.py database receives the full cohort"
        )

    def test_indel_database_output(self, module_text):
        body = _process_body(module_text, "RUN_DATABASE")
        assert "emit: indel_database" in body


# ===========================================================================
# 5. RUN_ANNOTATE
# ===========================================================================

class TestRunAnnotate:
    """indelible.py annotate requires --i, --o, --d, --config."""

    def test_annotate_command_flags(self, module_text):
        body = _process_body(module_text, "RUN_ANNOTATE")
        assert "indelible.py annotate" in body
        assert " --i " in body
        assert " --o " in body
        assert " --d " in body
        assert "--config" in body

    def test_annotated_output(self, module_text):
        body = _process_body(module_text, "RUN_ANNOTATE")
        assert "emit: annotated" in body


# ===========================================================================
# 6. INDELIBLE workflow – database broadcast via .first()
# ===========================================================================

class TestWorkflowDatabaseBroadcast:
    """The database channel must be converted to a value channel with .first()
    so that every scored file is annotated, not just the first one.

    Nextflow's default zip behaviour consumes one element from each input
    channel per invocation.  Because RUN_DATABASE emits exactly ONE database
    file while RUN_SCORE emits N scored files (one per sample), calling
    RUN_ANNOTATE(database_ch, scores_ch) without .first() would only annotate
    a single sample.  The .first() operator converts the single-element queue
    channel into a value channel that Nextflow reuses for every scores element.
    """

    def test_first_operator_used_on_database(self, module_text):
        """RUN_ANNOTATE must be called with .first() on the database channel."""
        wf_body = _workflow_body(module_text)
        assert ".first()" in wf_body, (
            "The INDELIBLE workflow must call RUN_ANNOTATE with "
            "RUN_DATABASE.out.indel_database.first() so the single database "
            "file is broadcast to every sample's scored file. Without .first(), "
            "only the first sample would be annotated."
        )

    def test_annotate_call_uses_first(self, module_text):
        """The RUN_ANNOTATE call must reference indel_database.first()."""
        wf_body = _workflow_body(module_text)
        assert re.search(r"indel_database\.first\(\)", wf_body), (
            "RUN_ANNOTATE must receive 'RUN_DATABASE.out.indel_database.first()' "
            "as its first argument"
        )


# ===========================================================================
# 7–9. De-novo processes
# ===========================================================================

class TestDenovoProcesses:
    """Each de-novo process must pass the correct parent BAM flags."""

    def test_trio_uses_both_parent_flags(self, module_text):
        """RUN_DENOVO_TRIO must pass both --m (maternal) and --p (paternal)."""
        body = _process_body(module_text, "RUN_DENOVO_TRIO")
        assert "indelible.py denovo" in body
        assert " --m " in body
        assert " --p " in body

    def test_mom_only_flag(self, module_text):
        """RUN_DENOVO_MOM must pass --m but not --p."""
        body = _process_body(module_text, "RUN_DENOVO_MOM")
        assert "indelible.py denovo" in body
        assert " --m " in body
        assert " --p " not in body, (
            "RUN_DENOVO_MOM should only pass --m (maternal); --p would require a "
            "paternal BAM that is not available in the mom-only denovo scenario."
        )

    def test_dad_only_flag(self, module_text):
        """RUN_DENOVO_DAD must pass --p but not --m."""
        body = _process_body(module_text, "RUN_DENOVO_DAD")
        assert "indelible.py denovo" in body
        assert " --p " in body
        assert " --m " not in body, (
            "RUN_DENOVO_DAD should only pass --p (paternal); --m would require a "
            "maternal BAM that is not available in the dad-only denovo scenario."
        )

    def test_all_denovo_use_annotation_flag(self, module_text):
        """All three denovo processes must pass --c (child annotation file)."""
        for proc in ("RUN_DENOVO_TRIO", "RUN_DENOVO_MOM", "RUN_DENOVO_DAD"):
            body = _process_body(module_text, proc)
            assert " --c " in body, f"{proc} must pass --c (child annotation) to indelible.py denovo"


# ===========================================================================
# 10. FILTER_INDELIBLE
# ===========================================================================

class TestFilterIndelibleProcess:
    """FILTER_INDELIBLE must apply the canonical awk AF/BP frequency filter."""

    def test_awk_filter_columns(self, module_text):
        """The awk filter must check columns $39 (AF_freq) and $40 (BP_freq)."""
        body = _process_body(module_text, "FILTER_INDELIBLE")
        assert "$39" in body and "$40" in body, (
            "FILTER_INDELIBLE awk must reference columns $39 (AF_freq) and $40 (BP_freq) "
            "from the 40-column annotated TSV produced by indelible.py annotate"
        )

    def test_sample_in_output_tuple(self, module_text):
        """The output must include val(sample) so downstream VCF conversion can
        use the sample ID to name output files correctly."""
        body = _process_body(module_text, "FILTER_INDELIBLE")
        output_section = re.search(r"output:\s*(.+?)(?=\s*script:)", body, re.DOTALL)
        assert output_section, "output section not found in FILTER_INDELIBLE"
        out_text = output_section.group(1)
        assert re.search(r"tuple\s+val\(sample\)", out_text), (
            "FILTER_INDELIBLE output must be 'tuple val(sample), path(...)' so that "
            "CONVERT_INDELIBLE_TO_VCF receives the sample ID for output file naming."
        )


# ===========================================================================
# 11–12. BGZIP_SORT_INDEX_VCF
# ===========================================================================

class TestBgzipSortIndexVcf:
    """BGZIP_SORT_INDEX_VCF must compress/sort/index only — no re-annotation."""

    def test_no_redundant_tool_annotation(self, module_text):
        """The process must NOT run a second bcftools-annotate for TOOL=INDELIBLE.

        CONVERT_INDELIBLE_TO_VCF already writes TOOL=INDELIBLE into every INFO
        field AND declares ##INFO=<ID=TOOL,...> in the VCF header.  Running
        bcftools annotate again with -h extra_header.txt would insert a duplicate
        ##INFO=<ID=TOOL,...> line, producing a non-standard VCF that downstream
        tools (SURVIVOR, Truvari) may reject or mishandle.
        """
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "bcftools annotate" not in body, (
            "BGZIP_SORT_INDEX_VCF must NOT run bcftools annotate to add TOOL=INDELIBLE "
            "because CONVERT_INDELIBLE_TO_VCF already writes that field. A second "
            "annotation would produce a duplicate ##INFO=<ID=TOOL,...> header line."
        )

    def test_no_extra_header_txt(self, module_text):
        """extra_header.txt must not be created; its purpose is already served by the converter."""
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "extra_header.txt" not in body, (
            "extra_header.txt adds a TOOL INFO header that already exists in the VCF "
            "from CONVERT_INDELIBLE_TO_VCF. Creating it again causes duplicate headers."
        )

    def test_bcftools_sort_present(self, module_text):
        """bcftools sort must be used to produce a coordinate-sorted compressed VCF."""
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "bcftools sort" in body

    def test_bgzip_output_flag(self, module_text):
        """bcftools sort must write bgzip-compressed output (-O z)."""
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "-O z" in body, (
            "bcftools sort must use '-O z' to write bgzip-compressed output directly"
        )

    def test_tabix_indexing(self, module_text):
        """tabix must be called to create the .tbi index for the sorted VCF."""
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "tabix" in body

    def test_sorted_vcf_gz_output(self, module_text):
        """Output must include *.sorted.vcf.gz and *.sorted.vcf.gz.tbi."""
        body = _process_body(module_text, "BGZIP_SORT_INDEX_VCF")
        assert "sorted.vcf.gz" in body
        assert "sorted.vcf.gz.tbi" in body


# ===========================================================================
# 13. INDELIBLE workflow – emit declarations
# ===========================================================================

class TestWorkflowEmits:
    """The INDELIBLE workflow must emit all expected channels."""

    @pytest.mark.parametrize("emit_name", [
        "filtered_cnvs",
        "indelible_denovo",
        "indelible_denovo_mom",
        "indelible_denovo_dad",
        "sorted_vcf",
        "sorted_vcf_index",
    ])
    def test_workflow_emits(self, module_text, emit_name):
        wf_body = _workflow_body(module_text)
        assert emit_name in wf_body, (
            f"INDELIBLE workflow must emit '{emit_name}' for downstream consumers"
        )
