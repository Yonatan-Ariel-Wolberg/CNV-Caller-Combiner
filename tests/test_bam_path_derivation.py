#!/usr/bin/env python3
"""
Tests for BAM-to-BAI path derivation and CNVkit BAM glob pattern in main.nf.

Validates two related workflow flaws that can cause pipeline failures:

  1. BAI path derivation for CANOES, XHMM and CLAMMS samplesheet workflows
     (main.nf case['canoes'], case['xhmm'], case['clamms']):
     The old regex "\\b.bam\\b" used an unescaped '.' (matches any character)
     and lacked a '$' end-of-string anchor, so a BAM path whose parent directory
     contained the substring '.bam.' (e.g. '/data/my.bam.dir/sample.bam') would
     have BOTH the directory component and the file extension replaced, producing
     a corrupt BAI path.
     The fix uses the Groovy slashy pattern /\\.bam$/ which matches only the
     literal '.bam' at the very end of the string.

  2. CNVkit BAM glob (main.nf case['cnvkit']):
     The old code used Java String.replace(".bam", ".{bam,bai}") which replaces
     ALL occurrences of the literal substring '.bam' in the path, not just the
     final extension.  A path like '/data/my.bam.dir/*.bam' would become
     '/data/my.{bam,bai}.dir/*.{bam,bai}', corrupting the directory component.
     The fix uses replaceAll(/\\.bam$/, '.{bam,bai}') which is end-anchored.
"""

import os
import re

import pytest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAIN_NF = os.path.join(REPO_ROOT, "main.nf")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_main():
    with open(MAIN_NF) as fh:
        return fh.read()


def _get_case_body(main_text, case_name):
    m = re.search(
        rf"case\['{re.escape(case_name)}'\](.+?)break",
        main_text,
        re.DOTALL,
    )
    assert m, f"case['{case_name}'] block not found in main.nf"
    return m.group(1)


# ---------------------------------------------------------------------------
# Python simulation of BAI-path derivation
# ---------------------------------------------------------------------------

def _derive_bai_path_fixed(bam_path):
    """Simulate the FIXED Groovy regex replaceAll(/\\.bam$/, '.bam.bai')."""
    return re.sub(r'\.bam$', '.bam.bai', bam_path)


def _derive_bai_path_buggy(bam_path):
    """Simulate the BUGGY Groovy regex replaceAll("\\b.bam\\b", ".bam.bai").

    In Java regex, '\\b.bam\\b' is: word-boundary + any-char + 'bam' + word-boundary.
    replaceAll() replaces ALL non-overlapping matches.
    """
    return re.sub(r'\b.bam\b', '.bam.bai', bam_path)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def main_text():
    return _read_main()


# ===========================================================================
# 1. Source-code checks: correct regex in main.nf
# ===========================================================================

class TestBaiRegexInSourceCode:
    """main.nf must use the end-anchored /\\.bam$/ pattern for BAI derivation."""

    def _check_case_uses_correct_bai_regex(self, main_text, case_name):
        body = _get_case_body(main_text, case_name)
        # Must NOT use the old buggy pattern (unescaped dot + word boundary).
        # In the .nf source file the old pattern is the literal text \\b.bam\\b
        # (with two backslashes), so we search for that exact byte sequence.
        assert '\\\\b.bam\\\\b' not in body, (
            f"case['{case_name}'] must not use '\\\\b.bam\\\\b' regex for BAI path "
            f"derivation; that pattern uses an unescaped '.' and lacks a '$' anchor"
        )
        # Must use an end-anchored replacement
        assert r'\.bam$' in body, (
            f"case['{case_name}'] must derive BAI path using an end-anchored regex "
            f"such as /\\.bam$/ so that only the final '.bam' extension is replaced"
        )

    def test_canoes_bai_regex_is_end_anchored(self, main_text):
        """case['canoes'] BAI path derivation must use /\\.bam$/ (not \\b.bam\\b)."""
        self._check_case_uses_correct_bai_regex(main_text, 'canoes')

    def test_xhmm_bai_regex_is_end_anchored(self, main_text):
        """case['xhmm'] BAI path derivation must use /\\.bam$/ (not \\b.bam\\b)."""
        self._check_case_uses_correct_bai_regex(main_text, 'xhmm')

    def test_clamms_bai_regex_is_end_anchored(self, main_text):
        """case['clamms'] BAI path derivation must use /\\.bam$/ (not \\b.bam\\b)."""
        self._check_case_uses_correct_bai_regex(main_text, 'clamms')

    def test_cnvkit_bam_glob_uses_end_anchored_replace(self, main_text):
        """case['cnvkit'] BAM glob must use replaceAll(/\\.bam$/) not replace('.bam')."""
        body = _get_case_body(main_text, 'cnvkit')
        # Must NOT use the old .replace() which replaces ALL occurrences
        assert '.replace(".bam"' not in body, (
            "case['cnvkit'] must not use .replace(\".bam\", ...) which replaces ALL "
            "occurrences of '.bam' in the path string, including directory components; "
            "use .replaceAll(/\\.bam$/, ...) to target only the file extension"
        )
        # Must use end-anchored replaceAll
        assert r'\.bam$' in body or ".replaceAll(" in body, (
            "case['cnvkit'] must use an end-anchored replaceAll to build the BAM/BAI "
            "glob pattern (e.g. .replaceAll(/\\.bam$/, '.{bam,bai}'))"
        )


# ===========================================================================
# 2. Behavioural correctness of the fixed BAI regex
# ===========================================================================

class TestBaiRegexBehaviour:
    """The fixed /\\.bam$/ regex must correctly derive BAI paths for all cases."""

    @pytest.mark.parametrize("bam_path,expected_bai", [
        # Standard paths
        ("/data/sample.bam",            "/data/sample.bam.bai"),
        ("/data/SAMPLE001.bam",         "/data/SAMPLE001.bam.bai"),
        ("sample.bam",                  "sample.bam.bai"),
        # Nested directories (no .bam in dir name)
        ("/home/user/project/s001.bam", "/home/user/project/s001.bam.bai"),
        # Sample ID contains 'bam' as a substring (not as extension)
        ("/data/mybam_sample.bam",      "/data/mybam_sample.bam.bai"),
    ])
    def test_fixed_regex_standard_paths(self, bam_path, expected_bai):
        """Fixed regex must correctly append .bai for common BAM path patterns."""
        assert _derive_bai_path_fixed(bam_path) == expected_bai, (
            f"Fixed BAI derivation failed for '{bam_path}'"
        )

    @pytest.mark.parametrize("bam_path,expected_bai", [
        # Path with '.bam.' in a directory component – the critical edge case
        ("/data/my.bam.dir/sample.bam",    "/data/my.bam.dir/sample.bam.bai"),
        ("/results/run.bam.archive/s.bam", "/results/run.bam.archive/s.bam.bai"),
    ])
    def test_fixed_regex_does_not_corrupt_directory_names(self, bam_path, expected_bai):
        """Fixed regex must NOT alter directory components that contain '.bam.'."""
        assert _derive_bai_path_fixed(bam_path) == expected_bai, (
            f"Fixed BAI derivation should leave directory components unchanged: "
            f"'{bam_path}' -> expected '{expected_bai}', "
            f"got '{_derive_bai_path_fixed(bam_path)}'"
        )

    @pytest.mark.parametrize("bam_path", [
        "/data/my.bam.dir/sample.bam",
        "/results/run.bam.archive/s.bam",
    ])
    def test_buggy_regex_corrupts_directory_names(self, bam_path):
        """Demonstrate that the OLD \\b.bam\\b regex corrupts paths with .bam.
        in directory names (proving the bug existed and why the fix is needed)."""
        buggy = _derive_bai_path_buggy(bam_path)
        fixed = _derive_bai_path_fixed(bam_path)
        # The buggy regex replaces MORE than just the final extension
        assert buggy != fixed, (
            f"Expected the buggy regex to produce a different (corrupt) result than "
            f"the fixed regex for '{bam_path}', but both returned '{fixed}'"
        )
        # The buggy result should contain an extra replacement in the directory part
        assert '.bam.bai.' in buggy or buggy.count('.bam.bai') > 1, (
            f"The buggy regex should have corrupted the directory component of "
            f"'{bam_path}', producing '{buggy}'"
        )


# ===========================================================================
# 3. Behavioural correctness of the fixed CNVkit glob pattern
# ===========================================================================

class TestCnvkitGlobBehaviour:
    """The fixed replaceAll(/\\.bam$/) must build correct BAM/BAI glob patterns."""

    def _make_cnvkit_glob_fixed(self, bams_param):
        """Simulate the FIXED Groovy: params.bams.replaceAll(/\\.bam$/, '.{bam,bai}')"""
        return re.sub(r'\.bam$', '.{bam,bai}', bams_param)

    def _make_cnvkit_glob_buggy(self, bams_param):
        """Simulate the BUGGY Groovy: params.bams.replace(".bam", ".{bam,bai}")"""
        return bams_param.replace(".bam", ".{bam,bai}")

    @pytest.mark.parametrize("bams_param,expected_glob", [
        ("/data/samples/*.bam",         "/data/samples/*.{bam,bai}"),
        ("/data/sample.bam",            "/data/sample.{bam,bai}"),
        ("/data/cohort/all_samples.bam", "/data/cohort/all_samples.{bam,bai}"),
    ])
    def test_fixed_glob_standard_paths(self, bams_param, expected_glob):
        """Fixed CNVkit glob must correctly transform standard BAM paths."""
        assert self._make_cnvkit_glob_fixed(bams_param) == expected_glob, (
            f"Fixed CNVkit glob failed for '{bams_param}'"
        )

    @pytest.mark.parametrize("bams_param,expected_glob", [
        # Path with '.bam.' in a directory component – the critical edge case
        ("/data/my.bam.dir/*.bam",      "/data/my.bam.dir/*.{bam,bai}"),
        ("/results/run.bam.arch/*.bam", "/results/run.bam.arch/*.{bam,bai}"),
    ])
    def test_fixed_glob_does_not_corrupt_directory_names(self, bams_param, expected_glob):
        """Fixed CNVkit glob must NOT modify directory components containing '.bam'."""
        assert self._make_cnvkit_glob_fixed(bams_param) == expected_glob, (
            f"Fixed CNVkit glob should leave directory components unchanged: "
            f"'{bams_param}' -> expected '{expected_glob}', "
            f"got '{self._make_cnvkit_glob_fixed(bams_param)}'"
        )

    @pytest.mark.parametrize("bams_param", [
        # These paths contain '.bam' literally in a directory component
        # (e.g., '.bam.dir') so .replace(".bam") corrupts the directory name
        "/data/my.bam.dir/*.bam",
        "/results/run.bam.archive/*.bam",
    ])
    def test_buggy_glob_corrupts_directory_names(self, bams_param):
        """Demonstrate that the OLD .replace() corrupts paths with '.bam' in
        directory names (proving the bug existed and why the fix is needed)."""
        buggy = self._make_cnvkit_glob_buggy(bams_param)
        fixed = self._make_cnvkit_glob_fixed(bams_param)
        # The buggy replace() modifies ALL occurrences, including directory components
        assert buggy != fixed, (
            f"Expected the buggy glob to produce a different (corrupt) result than "
            f"the fixed glob for '{bams_param}'"
        )
        # The buggy result should contain more than one {bam,bai} substitution
        assert buggy.count('.{bam,bai}') > 1, (
            f"The buggy glob should have replaced '.bam' in the directory component "
            f"of '{bams_param}', producing multiple replacements in '{buggy}'"
        )