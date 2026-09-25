#!/usr/bin/env python3
"""
Tests for modules/sv-mergers/modules-survivor.nf SURVIVOR merge parameters.

Verifies that the module contains exactly two SURVIVOR merge calls:
  1. Union merge        (min_support=1) → *_survivor_union.vcf
  2. Intersection merge (min_support=2) → *_survivor_intersection.vcf

Both calls share these fixed parameters:
  max_dist    = 1000  Maximum distance between breakpoints (bp)
  use_type    = 1     Require same SV type (1=yes)
  use_strand  = 0     CNV callers do not produce strand information
  est_dist    = 0     Do not estimate SV distance
  min_sv_size = 30    Minimum SV size (bp)
"""

import os
import re

MODULE_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'modules', 'sv-mergers', 'modules-survivor.nf'
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_module():
    with open(MODULE_PATH) as fh:
        return fh.read()


def _find_all_survivor_merge_args(content):
    """Return a list of parameter tuples for every SURVIVOR merge call found.

    In Nextflow script blocks backslashes are doubled (\\\\), so each
    line-continuation appears as \\\\ in the raw file string.
    Each tuple is (max_dist, min_support, use_type, use_strand, est_dist, min_sv_size).
    """
    pattern = re.compile(
        r'SURVIVOR\s+merge\s+\\\\\n'   # SURVIVOR merge \\
        r'\s+\S+\s+\\\\\n'             # \$list_file \\
        r'\s+(\d+)\s+\\\\\n'           # max_dist \\
        r'\s+(\d+)\s+\\\\\n'           # min_support \\
        r'\s+(\d+)\s+\\\\\n'           # use_type \\
        r'\s+(\d+)\s+\\\\\n'           # use_strand \\
        r'\s+(\d+)\s+\\\\\n'           # est_dist \\
        r'\s+(\d+)'                    # min_sv_size
    )
    return [tuple(int(g) for g in m.groups()) for m in pattern.finditer(content)]


# ---------------------------------------------------------------------------
# Tests – number of merge calls
# ---------------------------------------------------------------------------

def test_two_survivor_merge_calls():
    """Module must contain exactly two SURVIVOR merge calls (union + intersection)."""
    calls = _find_all_survivor_merge_args(_read_module())
    assert len(calls) == 2, (
        f"Expected 2 SURVIVOR merge calls (union + intersection), found {len(calls)}"
    )


# ---------------------------------------------------------------------------
# Tests – union merge (first call, min_support=1)
# ---------------------------------------------------------------------------

def test_union_merge_output_filename():
    """Union merge output filename must contain 'union'."""
    content = _read_module()
    assert '_survivor_union.vcf' in content, (
        "Expected union output filename '*_survivor_union.vcf' not found in module"
    )


def test_union_merge_min_support():
    """Union merge: min_support must be 1 (every SV from any caller is kept)."""
    calls = _find_all_survivor_merge_args(_read_module())
    _, min_support, *_ = calls[0]
    assert min_support == 1, f"Union min_support expected 1, got {min_support}"


# ---------------------------------------------------------------------------
# Tests – intersection merge (second call, min_support=2)
# ---------------------------------------------------------------------------

def test_intersection_merge_output_filename():
    """Intersection merge output filename must contain 'intersection'."""
    content = _read_module()
    assert '_survivor_intersection.vcf' in content, (
        "Expected intersection output filename '*_survivor_intersection.vcf' not found in module"
    )


def test_intersection_merge_min_support():
    """Intersection merge: min_support must be 2 (SV must be supported by >=2 callers)."""
    calls = _find_all_survivor_merge_args(_read_module())
    _, min_support, *_ = calls[1]
    assert min_support == 2, f"Intersection min_support expected 2, got {min_support}"


# ---------------------------------------------------------------------------
# Tests – shared parameters (same for both calls)
# ---------------------------------------------------------------------------

def test_both_calls_max_dist():
    """Both merges: max_dist must be 1000 bp."""
    for i, call in enumerate(_find_all_survivor_merge_args(_read_module())):
        max_dist = call[0]
        assert max_dist == 1000, f"Call {i}: max_dist expected 1000, got {max_dist}"


def test_both_calls_use_type():
    """Both merges: use_type must be 1 (require same SV type)."""
    for i, call in enumerate(_find_all_survivor_merge_args(_read_module())):
        use_type = call[2]
        assert use_type == 1, f"Call {i}: use_type expected 1, got {use_type}"


def test_both_calls_use_strand():
    """Both merges: use_strand must be 0 (CNV callers produce no strand information)."""
    for i, call in enumerate(_find_all_survivor_merge_args(_read_module())):
        use_strand = call[3]
        assert use_strand == 0, (
            f"Call {i}: use_strand expected 0 (CNV callers have no strand info), got {use_strand}"
        )


def test_both_calls_est_dist():
    """Both merges: est_dist must be 0."""
    for i, call in enumerate(_find_all_survivor_merge_args(_read_module())):
        est_dist = call[4]
        assert est_dist == 0, f"Call {i}: est_dist expected 0, got {est_dist}"


def test_both_calls_min_sv_size():
    """Both merges: min_sv_size must be 30 bp."""
    for i, call in enumerate(_find_all_survivor_merge_args(_read_module())):
        min_sv_size = call[5]
        assert min_sv_size == 30, f"Call {i}: min_sv_size expected 30, got {min_sv_size}"


# ---------------------------------------------------------------------------
# Tests – workflow emit names
# ---------------------------------------------------------------------------

def test_emit_union_vcf():
    """SURVIVOR workflow must emit 'union_vcf'."""
    content = _read_module()
    assert 'union_vcf' in content, "Expected 'union_vcf' emit in SURVIVOR workflow"


def test_emit_intersection_vcf():
    """SURVIVOR workflow must emit 'intersection_vcf'."""
    content = _read_module()
    assert 'intersection_vcf' in content, "Expected 'intersection_vcf' emit in SURVIVOR workflow"

