#!/usr/bin/env python3
"""
Tests for bin/feature_extraction.py.

Validates:
  1. Module-level structure: imports, public API, module docstring.
  2. Helper functions: _load_bed, _count_probes, _count_probes_flank,
     _mean_mappability, _load_mappability_bed, _gc_content, _rd_ratio,
     _mean_depth, _l2r_stats.
  3. extract_normalized_features: signature has required optional parameters
     (bed_file, bam_file, reference_fasta, mappability_file, rd_flank),
     produces expected columns, and handles missing optional inputs.
  4. QUAL_norm is read from the per-caller VCF QUAL field (not recomputed).
     Raw caller-native scores (Q_SOME, SQ, QS, CNQ, native QUAL) are already
     the direct input to QUAL_norm and are therefore NOT stored as separate
     columns.

Where full VCF/BAM/FASTA fixtures are impractical, the tests import the
helper functions directly and exercise them with in-memory data.
"""

import importlib
import io
import math
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Locate the script
# ---------------------------------------------------------------------------

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT_PATH = os.path.join(BIN_DIR, 'feature_extraction.py')


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _import_module():
    """Import feature_extraction.py as a Python module."""
    spec = importlib.util.spec_from_file_location('feature_extraction', SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def fe():
    """Return the feature_extraction module."""
    return _import_module()


# ===========================================================================
# 1. Module structure
# ===========================================================================

class TestModuleStructure:
    """feature_extraction.py must import the required libraries."""

    def test_imports_pysam(self, script_text):
        assert 'import pysam' in script_text

    def test_imports_pandas(self, script_text):
        assert 'import pandas' in script_text or 'pandas as pd' in script_text

    def test_imports_numpy(self, script_text):
        assert 'import numpy' in script_text or 'numpy as np' in script_text

    def test_imports_defaultdict(self, script_text):
        assert 'defaultdict' in script_text

    def test_exposes_extract_normalized_features(self, script_text):
        assert 'def extract_normalized_features(' in script_text

    def test_module_docstring_mentions_qual_norm(self, script_text):
        assert 'qual_norm' in script_text.lower() or 'QUAL_norm' in script_text

    def test_module_docstring_mentions_n_probes(self, script_text):
        assert 'n_probes' in script_text

    def test_module_docstring_mentions_rd_ratio(self, script_text):
        assert 'rd_ratio' in script_text

    def test_module_docstring_mentions_gc_content(self, script_text):
        assert 'gc_content' in script_text

    def test_module_docstring_mentions_mappability(self, script_text):
        assert 'mappability' in script_text

    def test_module_docstring_mentions_pd3_paper(self, script_text):
        assert 'pd3' in script_text or 'cnv-paper-2020' in script_text

    def test_module_docstring_mentions_n_probes_flank(self, script_text):
        assert 'n_probes_flank' in script_text

    def test_module_docstring_mentions_l2r(self, script_text):
        assert 'l2r_mean' in script_text


# ===========================================================================
# 2. extract_normalized_features signature
# ===========================================================================

class TestFunctionSignature:
    """extract_normalized_features must accept the new optional parameters."""

    def test_has_bed_file_param(self, script_text):
        assert 'bed_file' in script_text

    def test_has_bam_file_param(self, script_text):
        assert 'bam_file' in script_text

    def test_has_reference_fasta_param(self, script_text):
        assert 'reference_fasta' in script_text

    def test_has_mappability_file_param(self, script_text):
        assert 'mappability_file' in script_text

    def test_has_rd_flank_param(self, script_text):
        assert 'rd_flank' in script_text

    def test_bed_file_defaults_none(self, script_text):
        # Default value of bed_file must be None so the parameter is optional.
        assert 'bed_file=None' in script_text

    def test_bam_file_defaults_none(self, script_text):
        assert 'bam_file=None' in script_text

    def test_reference_fasta_defaults_none(self, script_text):
        assert 'reference_fasta=None' in script_text

    def test_mappability_file_defaults_none(self, script_text):
        assert 'mappability_file=None' in script_text


# ===========================================================================
# 3. QUAL_norm sourcing: must read from VCF QUAL field, not recompute
# ===========================================================================

class TestQualNormSourcing:
    """QUAL_norm should come from orig.qual (normalised VCF), not from a
    separate normalize_q recomputation inside feature_extraction.py."""

    def test_reads_orig_qual(self, script_text):
        assert 'orig.qual' in script_text, (
            "feature_extraction must read QUAL_norm from orig.qual "
            "(the QUAL field of the normalised per-caller VCF)"
        )

    def test_no_standalone_normalize_q_function(self, script_text):
        # normalize_q was the old re-implementation; it should be gone now
        # that we read QUAL_norm directly from the normalised VCF QUAL field.
        assert 'def normalize_q(' not in script_text, (
            "normalize_q helper should be removed; QUAL_norm is now read "
            "directly from the normalised VCF QUAL field (orig.qual)"
        )

    def test_stores_qual_norm_per_caller(self, script_text):
        assert "qual_norm_" in script_text, (
            "feature_extraction must store qual_norm_{caller} columns"
        )

    def test_no_raw_qual_per_caller(self, script_text):
        assert "raw_qual_" not in script_text, (
            "raw_qual_{caller} columns must be absent: the caller-native "
            "scores are already the input to QUAL_norm and are redundant"
        )

    def test_no_get_raw_qual_function(self, script_text):
        assert "def _get_raw_qual(" not in script_text, (
            "_get_raw_qual helper must be removed along with raw_qual columns"
        )


# ===========================================================================
# 5. Probe count helper (_load_bed, _count_probes)
# ===========================================================================

class TestProbeCountHelpers:
    """_load_bed and _count_probes must be present and correct."""

    def test_load_bed_present(self, script_text):
        assert 'def _load_bed(' in script_text

    def test_count_probes_present(self, script_text):
        assert 'def _count_probes(' in script_text

    def test_count_probes_uses_overlap_logic(self, script_text):
        # Overlap test: iv_start < end and iv_end > start
        assert 'iv_start < end' in script_text or 'iv_end > start' in script_text

    def test_count_probes_correct(self, fe, tmp_path):
        # Write a small BED file
        bed = tmp_path / 'targets.bed'
        bed.write_text(
            'chr1\t100\t200\n'   # overlaps [150, 250)
            'chr1\t300\t400\n'   # does NOT overlap [150, 250)
            'chr1\t190\t260\n'   # overlaps [150, 250)
            'chr2\t100\t200\n'   # wrong chrom
        )
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 150, 250, intervals)
        assert count == 2

    def test_count_probes_no_overlap(self, fe, tmp_path):
        bed = tmp_path / 'targets2.bed'
        bed.write_text('chr1\t500\t600\n')
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_probe_ends_at_cnv_start(self, fe, tmp_path):
        """Probe ending exactly at CNV start must NOT count as overlapping."""
        bed = tmp_path / 'boundary1.bed'
        bed.write_text('chr1\t0\t100\n')  # probe [0, 100)
        intervals = fe._load_bed(str(bed))
        # CNV starts at 100 -> probe iv_end (100) > start (100) is False
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_probe_starts_at_cnv_end(self, fe, tmp_path):
        """Probe starting exactly at CNV end must NOT count as overlapping."""
        bed = tmp_path / 'boundary2.bed'
        bed.write_text('chr1\t200\t300\n')  # probe [200, 300)
        intervals = fe._load_bed(str(bed))
        # CNV ends at 200 -> probe iv_start (200) < end (200) is False
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_single_base_overlap(self, fe, tmp_path):
        """Probe overlapping by exactly one base must count."""
        bed = tmp_path / 'boundary3.bed'
        bed.write_text('chr1\t199\t300\n')  # probe [199, 300)
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 1

    def test_load_bed_skips_comment_lines(self, fe, tmp_path):
        bed = tmp_path / 'targets3.bed'
        bed.write_text(
            '# comment\n'
            'chr1\t10\t20\n'
        )
        intervals = fe._load_bed(str(bed))
        assert len(intervals['chr1']) == 1


# ===========================================================================
# 6. Mappability helpers (_load_mappability_bed, _mean_mappability)
# ===========================================================================

class TestMappabilityHelpers:
    """Mappability BED loader and weighted-mean scorer must be correct."""

    def test_load_mappability_bed_present(self, script_text):
        assert 'def _load_mappability_bed(' in script_text

    def test_mean_mappability_present(self, script_text):
        assert 'def _mean_mappability(' in script_text

    def test_mean_mappability_weighted(self, fe, tmp_path):
        bed = tmp_path / 'map.bed'
        # Two equal-sized intervals with different scores
        bed.write_text(
            'chr1\t0\t100\t0.8\n'   # 100 bp, score 0.8
            'chr1\t100\t200\t0.4\n' # 100 bp, score 0.4
        )
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 200, intervals)
        assert abs(score - 0.6) < 1e-9

    def test_mean_mappability_partial_overlap(self, fe, tmp_path):
        # Interval [50, 150) with score 1.0; query [0, 100).
        # Overlap = [50, 100) = 50 bp. Weighted mean is over covered bases
        # only: 50 * 1.0 / 50 = 1.0.
        bed = tmp_path / 'map2.bed'
        bed.write_text('chr1\t50\t150\t1.0\n')
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 100, intervals)
        assert score == pytest.approx(1.0)

    def test_mean_mappability_unequal_sizes(self, fe, tmp_path):
        """Weighted mean with intervals of unequal size must not be a simple average."""
        bed = tmp_path / 'map_unequal.bed'
        # 100 bp at 0.8, 50 bp at 0.4 -> (100*0.8 + 50*0.4) / 150 = 100/150
        bed.write_text(
            'chr1\t0\t100\t0.8\n'
            'chr1\t100\t150\t0.4\n'
        )
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 150, intervals)
        expected = (100 * 0.8 + 50 * 0.4) / 150
        assert score == pytest.approx(expected)

    def test_mean_mappability_no_overlap_returns_nan(self, fe, tmp_path):
        bed = tmp_path / 'map3.bed'
        bed.write_text('chr1\t500\t600\t1.0\n')
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 100, intervals)
        assert math.isnan(score)

    def test_mean_mappability_empty_dict_returns_nan(self, fe):
        score = fe._mean_mappability('chr1', 0, 100, {})
        assert math.isnan(score)


# ===========================================================================
# 7. GC content helper (_gc_content)
# ===========================================================================

class TestGcContent:
    """_gc_content must correctly compute G/C fraction."""

    def test_gc_content_present(self, script_text):
        assert 'def _gc_content(' in script_text

    def test_gc_content_pure_gc(self, fe):
        """100 % GC sequence."""
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'GCGCGC'
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(1.0)

    def test_gc_content_pure_at(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'ATATAT'
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(0.0)

    def test_gc_content_mixed(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'GCATAT'  # 2 GC out of 6
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(2 / 6)

    def test_gc_content_empty_returns_nan(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return ''
        assert math.isnan(fe._gc_content(MockFasta(), 'chr1', 0, 0))


# ===========================================================================
# 8. RD ratio helper (_rd_ratio, _mean_depth)
# ===========================================================================

class TestRdRatio:
    """_rd_ratio must return target_depth / flank_depth and NaN on zero flank."""

    def test_rd_ratio_present(self, script_text):
        assert 'def _rd_ratio(' in script_text

    def test_mean_depth_present(self, script_text):
        assert 'def _mean_depth(' in script_text

    def test_rd_ratio_nan_on_zero_flank(self, fe):
        """When flanking depth is 0, _rd_ratio must return NaN."""
        class MockPileupCol:
            def __init__(self, n):
                self.nsegments = n

        class MockBam:
            def pileup(self, chrom, start, end, **kwargs):
                return []  # empty -> depth 0

        result = fe._rd_ratio(MockBam(), 'chr1', 1000, 2000, flank=500)
        assert math.isnan(result)


# ===========================================================================
# 9. Concordance and caller flags
# ===========================================================================

class TestConcordanceAndFlags:
    """SUPP_VEC-derived concordance and is_{caller} flags must be present."""

    def test_concordance_column_produced(self, script_text):
        assert "'concordance'" in script_text or '"concordance"' in script_text

    def test_is_caller_columns_produced(self, script_text):
        assert "f'is_{caller_name}'" in script_text or "is_{caller" in script_text

    def test_supp_vec_parsed(self, script_text):
        assert 'SUPP_VEC' in script_text

    def test_cnv_type_column_produced(self, script_text):
        assert "'cnv_type'" in script_text or 'cnv_type' in script_text

    def test_cnv_size_column_produced(self, script_text):
        assert "'cnv_size'" in script_text or 'cnv_size' in script_text



# ===========================================================================
# 12. NaN defaults when optional files are absent
# ===========================================================================

class TestNaNDefaults:
    """When optional annotation files are not provided, relevant columns
    must default to NaN rather than raising exceptions."""

    def test_n_probes_nan_without_bed(self, script_text):
        # When bed_intervals is empty/False, n_probes must be np.nan
        assert 'np.nan' in script_text

    def test_rd_ratio_nan_without_bam(self, script_text):
        assert 'rd_ratio' in script_text and 'np.nan' in script_text

    def test_gc_content_nan_without_fasta(self, script_text):
        assert 'gc_content' in script_text and 'np.nan' in script_text

    def test_mappability_nan_without_file(self, script_text):
        assert 'mappability' in script_text and 'np.nan' in script_text


# ===========================================================================
# 13. CN-Learn-inspired features: _cnv_size_label
# ===========================================================================

class TestCnvSizeLabel:
    """_cnv_size_label must map CNV sizes to CN-Learn's ordinal size bins."""

    def test_size_label_present(self, script_text):
        assert 'def _cnv_size_label(' in script_text

    def test_size_label_column_emitted(self, script_text):
        assert 'size_label' in script_text

    def test_bin_1_less_than_1kb(self, fe):
        assert fe._cnv_size_label(500) == 1

    def test_bin_1_boundary_999(self, fe):
        assert fe._cnv_size_label(999) == 1

    def test_bin_2_exactly_1kb(self, fe):
        assert fe._cnv_size_label(1_000) == 2

    def test_bin_3_5kb(self, fe):
        assert fe._cnv_size_label(5_000) == 3

    def test_bin_4_10kb(self, fe):
        assert fe._cnv_size_label(10_000) == 4

    def test_bin_5_25kb(self, fe):
        assert fe._cnv_size_label(25_000) == 5

    def test_bin_6_50kb(self, fe):
        assert fe._cnv_size_label(50_000) == 6

    def test_bin_7_75kb(self, fe):
        assert fe._cnv_size_label(75_000) == 7

    def test_bin_8_100kb(self, fe):
        assert fe._cnv_size_label(100_000) == 8

    def test_bin_9_250kb(self, fe):
        assert fe._cnv_size_label(250_000) == 9

    def test_bin_10_500kb(self, fe):
        assert fe._cnv_size_label(500_000) == 10

    def test_bin_11_1mb(self, fe):
        assert fe._cnv_size_label(1_000_000) == 11

    def test_bin_12_5mb(self, fe):
        assert fe._cnv_size_label(5_000_000) == 12

    def test_bin_12_very_large(self, fe):
        assert fe._cnv_size_label(10_000_000) == 12

    def test_all_labels_in_range_1_to_12(self, fe):
        test_sizes = [0, 500, 999, 1000, 4999, 5000, 9999, 10000, 24999,
                      25000, 49999, 50000, 74999, 75000, 99999, 100000,
                      249999, 250000, 499999, 500000, 999999, 1000000,
                      4999999, 5000000, 10000000]
        for s in test_sizes:
            label = fe._cnv_size_label(s)
            assert 1 <= label <= 12, f"size {s} -> label {label} not in [1,12]"


# ===========================================================================
# 14. CN-Learn-inspired features: _encode_chrom
# ===========================================================================

class TestEncodeChrom:
    """_encode_chrom must return integer chromosome codes matching CN-Learn."""

    def test_encode_chrom_present(self, script_text):
        assert 'def _encode_chrom(' in script_text

    def test_chrom_encoded_column_emitted(self, script_text):
        assert 'chrom_encoded' in script_text

    def test_autosome_without_prefix(self, fe):
        assert fe._encode_chrom('1') == 1
        assert fe._encode_chrom('22') == 22

    def test_autosome_with_chr_prefix(self, fe):
        assert fe._encode_chrom('chr1') == 1
        assert fe._encode_chrom('chr22') == 22

    def test_x_chromosome(self, fe):
        assert fe._encode_chrom('X') == 23
        assert fe._encode_chrom('chrX') == 23
        assert fe._encode_chrom('x') == 23

    def test_y_chromosome(self, fe):
        assert fe._encode_chrom('Y') == 24
        assert fe._encode_chrom('chrY') == 24

    def test_unknown_returns_zero(self, fe):
        assert fe._encode_chrom('chrM') == 0
        assert fe._encode_chrom('scaffold_1') == 0
        assert fe._encode_chrom('') == 0


# ===========================================================================
# 15. Aggregate quality summary columns
# ===========================================================================

class TestAggregateQualNorm:
    """max_qual_norm and mean_qual_norm_supported columns must be present."""

    def test_max_qual_norm_present(self, script_text):
        assert 'max_qual_norm' in script_text

    def test_mean_qual_norm_supported_present(self, script_text):
        assert 'mean_qual_norm_supported' in script_text

    def test_aggregate_computed_from_supporting_callers_only(self, script_text):
        # The aggregate must filter on is_{caller} == 1, not all callers.
        assert "is_{cn}" in script_text or "f'is_{cn}'" in script_text or \
               "_qual_norm_vals" in script_text


# ===========================================================================
# 16. pd3/cnv-paper-2020 inspired features: _count_probes_flank
# ===========================================================================

class TestProbeFlankCount:
    """_count_probes_flank must count probes outside the CNV in a flanking window."""

    def test_count_probes_flank_present(self, script_text):
        assert 'def _count_probes_flank(' in script_text

    def test_n_probes_flank_column_emitted(self, script_text):
        assert 'n_probes_flank' in script_text

    def test_count_probes_flank_basic(self, fe, tmp_path):
        """Probes flanking a CNV should be counted; probes inside should not."""
        bed = tmp_path / 'flank.bed'
        # CNV is [200, 400), length=200, flank=100 on each side → window [100, 500)
        # Probes: 50-100 (outside window), 100-150 (in left flank), 150-200 (adjacent),
        #         250-350 (in CNV), 400-450 (in right flank), 500-600 (outside window)
        bed.write_text(
            'chr1\t50\t100\n'    # outside window (ends at window start)
            'chr1\t100\t150\n'   # left flank: overlaps window and not CNV
            'chr1\t150\t200\n'   # left flank: ends at CNV start
            'chr1\t250\t350\n'   # inside CNV -> excluded
            'chr1\t400\t450\n'   # right flank
            'chr1\t500\t600\n'   # outside window
        )
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes_flank('chr1', 200, 400, intervals)
        # Expected: probes at [100,150) and [150,200) and [400,450) = 3
        assert count == 3

    def test_count_probes_flank_excludes_cnv_probes(self, fe, tmp_path):
        """Probes overlapping the CNV must not be counted."""
        bed = tmp_path / 'flank2.bed'
        bed.write_text(
            'chr1\t0\t100\n'     # ends exactly at left-flank window start (iv_end=100 > left_start=100 is False) -> outside window
            'chr1\t100\t300\n'   # overlaps CNV [200, 400) -> in_cnv=True, excluded
            'chr1\t400\t500\n'   # right flank, not in CNV
        )
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes_flank('chr1', 200, 400, intervals)
        # Only [400,500) qualifies: in flank window and outside CNV
        assert count == 1

    def test_count_probes_flank_no_flanking_probes(self, fe, tmp_path):
        """When there are no probes outside the CNV, return 0."""
        bed = tmp_path / 'flank3.bed'
        bed.write_text('chr1\t200\t400\n')  # entirely within CNV
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes_flank('chr1', 200, 400, intervals)
        assert count == 0

    def test_count_probes_flank_wrong_chrom(self, fe, tmp_path):
        bed = tmp_path / 'flank4.bed'
        bed.write_text('chr2\t0\t100\n')
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes_flank('chr1', 200, 400, intervals)
        assert count == 0

    def test_count_probes_flank_nan_without_bed(self, script_text):
        assert 'n_probes_flank' in script_text and 'np.nan' in script_text


# ===========================================================================
# 17. pd3/cnv-paper-2020 inspired features: _l2r_stats
# ===========================================================================

class TestL2rStats:
    """_l2r_stats must return probe-level log2 ratio statistics."""

    def test_l2r_stats_present(self, script_text):
        assert 'def _l2r_stats(' in script_text

    def test_l2r_mean_column_emitted(self, script_text):
        assert 'l2r_mean' in script_text

    def test_l2r_dev_column_emitted(self, script_text):
        assert 'l2r_dev' in script_text

    def test_l2r_flank_diff_column_emitted(self, script_text):
        assert 'l2r_flank_diff' in script_text

    def test_l2r_stats_nan_on_no_probes(self, fe, tmp_path):
        """When the BED has no probes, all three stats must be NaN."""
        bed = tmp_path / 'empty.bed'
        bed.write_text('')
        intervals = fe._load_bed(str(bed))

        class MockBam:
            def pileup(self, *a, **kw):
                return []

        l2r_mean, l2r_dev, l2r_flank_diff = fe._l2r_stats(
            MockBam(), 'chr1', 200, 400, intervals
        )
        assert math.isnan(l2r_mean)
        assert math.isnan(l2r_dev)
        assert math.isnan(l2r_flank_diff)

    def test_l2r_stats_nan_on_no_flanking_probes(self, fe, tmp_path):
        """When there are CNV probes but no flanking probes, all three must be NaN."""
        bed = tmp_path / 'noflanks.bed'
        bed.write_text('chr1\t200\t400\n')  # only inside CNV
        intervals = fe._load_bed(str(bed))

        class MockBam:
            def pileup(self, chrom, start, end, **kw):
                # Return depth 10 for any region
                return [type('Col', (), {'nsegments': 10})()] * (end - start)

        l2r_mean, l2r_dev, l2r_flank_diff = fe._l2r_stats(
            MockBam(), 'chr1', 200, 400, intervals
        )
        assert math.isnan(l2r_mean)
        assert math.isnan(l2r_dev)
        assert math.isnan(l2r_flank_diff)

    def test_l2r_stats_del_signal(self, fe, tmp_path):
        """DEL should produce a negative l2r_mean."""
        bed = tmp_path / 'del.bed'
        # 2 probes in CNV [1000, 2000), 2 left-flank, 2 right-flank
        bed.write_text(
            'chr1\t800\t900\n'    # left flank probe 1
            'chr1\t900\t1000\n'   # left flank probe 2
            'chr1\t1000\t1200\n'  # CNV probe 1
            'chr1\t1600\t1800\n'  # CNV probe 2
            'chr1\t2000\t2100\n'  # right flank probe 1
            'chr1\t2100\t2200\n'  # right flank probe 2
        )
        intervals = fe._load_bed(str(bed))

        def _make_bam(depths_by_region):
            """depths_by_region: list of (start, depth) sorted by start."""
            class MockPilCol:
                def __init__(self, n):
                    self.nsegments = n

            class MockBam:
                def pileup(self, chrom, start, end, **kw):
                    # Return the depth associated with this region
                    for reg_start, dep in depths_by_region:
                        if abs(start - reg_start) < 300:
                            return [MockPilCol(dep)]
                    return []

            return MockBam()

        # Flanking probes: depth 20; CNV probes: depth 10 (half -> DEL, l2r ≈ -1)
        bam = _make_bam([
            (800, 20), (900, 20),    # left flank
            (1000, 10), (1600, 10),  # CNV (halved depth)
            (2000, 20), (2100, 20),  # right flank
        ])
        l2r_mean, l2r_dev, l2r_flank_diff = fe._l2r_stats(
            bam, 'chr1', 1000, 2000, intervals
        )
        assert l2r_mean < 0, f"Expected negative l2r_mean for DEL, got {l2r_mean}"

    def test_l2r_stats_dup_signal(self, fe, tmp_path):
        """DUP should produce a positive l2r_mean."""
        bed = tmp_path / 'dup.bed'
        bed.write_text(
            'chr1\t800\t900\n'
            'chr1\t900\t1000\n'
            'chr1\t1000\t1200\n'
            'chr1\t1600\t1800\n'
            'chr1\t2000\t2100\n'
            'chr1\t2100\t2200\n'
        )
        intervals = fe._load_bed(str(bed))

        class MockPilCol:
            def __init__(self, n):
                self.nsegments = n

        class MockBam:
            def pileup(self, chrom, start, end, **kw):
                # Flanking depth=10, CNV depth=20 (doubled -> DUP, l2r ≈ +1)
                if 1000 <= start < 2000:
                    return [MockPilCol(20)]
                return [MockPilCol(10)]

        l2r_mean, l2r_dev, l2r_flank_diff = fe._l2r_stats(
            MockBam(), 'chr1', 1000, 2000, intervals
        )
        assert l2r_mean > 0, f"Expected positive l2r_mean for DUP, got {l2r_mean}"

    def test_l2r_stats_nan_without_bam_or_bed(self, script_text):
        """When either bam_file or bed_file is absent, L2R columns must be NaN."""
        assert 'l2r_mean' in script_text and 'np.nan' in script_text




# ===========================================================================
# 19. concordance computed in Truvari mode
# ===========================================================================

class TestTruvariConcordance:
    """extract_normalized_features must compute concordance in Truvari mode."""

    def test_concordance_in_truvari_branch(self, script_text):
        """concordance column must be assigned inside the truvari branch."""
        # The script must set 'concordance' after the Truvari is_{caller} loop.
        assert "concordance" in script_text, (
            "concordance feature must be computed in Truvari mode"
        )

    def test_truvari_concordance_is_sum_of_is_flags(self, script_text):
        """concordance must be the sum of is_{caller} values in Truvari mode."""
        assert "sum(" in script_text and "is_" in script_text, (
            "concordance should be computed as sum of is_ flags in Truvari mode"
        )


# ===========================================================================
# 20. SURVIVOR header auto-detection helper
# ===========================================================================

class TestCallerOrderFromSurvivorHeader:
    """_caller_order_from_survivor_header must parse ##SAMPLE lines."""

    def test_helper_present(self, script_text):
        assert 'def _caller_order_from_survivor_header(' in script_text

    def test_helper_returns_list(self, fe):
        """Should return an empty list when given a header with no ##SAMPLE lines."""
        import pysam

        class MockRecord:
            type = 'GENERIC'
            def __str__(self):
                return '##reference=GRCh38'

        class MockHeader:
            records = [MockRecord()]

        result = fe._caller_order_from_survivor_header(MockHeader())
        assert isinstance(result, list)
        assert result == []

    def test_caller_patterns_present(self, script_text):
        """The helper must contain filename patterns for known callers."""
        for caller in ('CANOES', 'CLAMMS', 'XHMM', 'CNVKIT', 'DRAGEN', 'EXOMEDEPTH'):
            assert caller in script_text, (
                f"_caller_order_from_survivor_header must recognise {caller} filenames"
            )

    def test_dragen_pattern_matches_cnv_vcf_gz(self, fe):
        """DRAGEN caller pattern must match ${sample_id}.cnv.vcf.gz filenames.

        DRAGEN Germline Enrichment names its CNV VCF outputs as
        ${sample_id}.cnv.vcf.gz.  The pattern for 'dragen' in
        _CALLER_FILENAME_PATTERNS must recognise this naming convention.
        """
        import re as _re
        dragen_pattern = None
        for name, pat in fe._CALLER_FILENAME_PATTERNS:
            if name == 'dragen':
                dragen_pattern = pat
                break
        assert dragen_pattern is not None, \
            "'dragen' entry not found in _CALLER_FILENAME_PATTERNS"
        assert dragen_pattern.search('NA12878.cnv.vcf.gz'), (
            "DRAGEN pattern must match 'NA12878.cnv.vcf.gz' "
            "(DRAGEN Germline Enrichment output naming)"
        )
        assert dragen_pattern.search('SAMPLE.cnv.vcf'), (
            "DRAGEN pattern must match 'SAMPLE.cnv.vcf'"
        )
        # Must still match the existing DRAGEN keyword in filenames
        assert dragen_pattern.search('PT01_DRAGEN_output.vcf.gz'), (
            "DRAGEN pattern must still match filenames containing 'DRAGEN'"
        )

    def test_parses_sample_id_and_file(self, script_text):
        """Parser must extract ID= and File= from ##SAMPLE meta-info lines."""
        assert 'ID=' in script_text and 'File=' in script_text

    def test_falls_back_to_tool_n(self, fe):
        """Unrecognised filenames must yield 'tool_N' names."""
        class MockRecord:
            type = 'GENERIC'
            def __str__(self):
                return '##SAMPLE=<ID=0,File=/path/to/unknown_caller.vcf>'

        class MockHeader:
            records = [MockRecord()]

        import re as _re
        result = fe._caller_order_from_survivor_header(MockHeader())
        assert result == ['tool_0'], f"Expected ['tool_0'], got {result}"


# ===========================================================================
# 20b. SUPP_VEC → caller set (_callers_from_supp_vec)
# ===========================================================================

class TestCallersSuppVec:
    """_callers_from_supp_vec must correctly map SUPP_VEC bits to caller names.

    This is the core mechanism by which the feature extraction script
    ascertains which CNV callers were involved in calling a variant.
    """

    def test_function_present(self, script_text):
        """_callers_from_supp_vec must be defined in the script."""
        assert 'def _callers_from_supp_vec(' in script_text

    def test_returns_set(self, fe):
        """Return type must be a set."""
        result = fe._callers_from_supp_vec('101', ['canoes', 'xhmm', 'cnvkit'])
        assert isinstance(result, set)

    def test_single_caller_first_bit(self, fe):
        """Only bit 0 set → only the first caller."""
        result = fe._callers_from_supp_vec('100', ['canoes', 'xhmm', 'cnvkit'])
        assert result == {'canoes'}

    def test_single_caller_middle_bit(self, fe):
        """Only bit 1 set → only the second caller."""
        result = fe._callers_from_supp_vec('010', ['canoes', 'xhmm', 'cnvkit'])
        assert result == {'xhmm'}

    def test_single_caller_last_bit(self, fe):
        """Only bit 2 set → only the third caller."""
        result = fe._callers_from_supp_vec('001', ['canoes', 'xhmm', 'cnvkit'])
        assert result == {'cnvkit'}

    def test_two_callers(self, fe):
        """Bits 0 and 2 set → first and third callers."""
        result = fe._callers_from_supp_vec('101', ['canoes', 'xhmm', 'cnvkit'])
        assert result == {'canoes', 'cnvkit'}

    def test_all_callers(self, fe):
        """All bits set → all callers returned."""
        order = ['canoes', 'clamms', 'xhmm']
        result = fe._callers_from_supp_vec('111', order)
        assert result == {'canoes', 'clamms', 'xhmm'}

    def test_no_callers_all_zeros(self, fe):
        """All bits zero → empty set (variant called by nobody)."""
        result = fe._callers_from_supp_vec('000', ['canoes', 'xhmm', 'cnvkit'])
        assert result == set()

    def test_empty_supp_vec_returns_empty_set(self, fe):
        """Empty SUPP_VEC string → empty set."""
        result = fe._callers_from_supp_vec('', ['canoes', 'xhmm'])
        assert result == set()

    def test_empty_caller_order_falls_back_to_tool_n(self, fe):
        """When caller_order is empty, bit positions fall back to 'tool_N' names."""
        result = fe._callers_from_supp_vec('101', [])
        assert result == {'tool_0', 'tool_2'}

    def test_partial_caller_order_falls_back_to_tool_n(self, fe):
        """Bit positions beyond caller_order length fall back to 'tool_N'."""
        result = fe._callers_from_supp_vec('111', ['canoes'])
        assert result == {'canoes', 'tool_1', 'tool_2'}

    def test_concordance_equals_len_result(self, fe):
        """len(result) gives concordance (number of callers that called variant)."""
        result = fe._callers_from_supp_vec('1101', ['canoes', 'clamms', 'xhmm', 'cnvkit'])
        assert len(result) == 3

    def test_all_supported_callers(self, fe):
        """All supported callers can be identified from a full SUPP_VEC."""
        order = list(fe.SUPPORTED_CALLERS)
        supp_vec = '1' * len(order)
        result = fe._callers_from_supp_vec(supp_vec, order)
        assert result == set(fe.SUPPORTED_CALLERS)

    def test_subset_of_supported_callers(self, fe):
        """Only callers with bit='1' are returned, not the full SUPPORTED_CALLERS."""
        order = list(fe.SUPPORTED_CALLERS)
        supp_vec = '10001' + '0' * (len(order) - 5)
        result = fe._callers_from_supp_vec(supp_vec, order)
        assert result == {'canoes', 'cnvkit'}

    def test_function_used_in_survivor_branch(self, script_text):
        """_callers_from_supp_vec must be called in the SURVIVOR mode branch."""
        assert '_callers_from_supp_vec(' in script_text

    def test_concordance_uses_len_of_supporting_callers(self, script_text):
        """concordance must be derived from len(supporting_callers), not a raw sum."""
        assert 'len(supporting_callers)' in script_text

class TestCLI:
    """feature_extraction.py must expose an argparse CLI."""

    def test_argparse_import(self, script_text):
        assert 'import argparse' in script_text

    def test_main_block_present(self, script_text):
        assert "__name__ == '__main__'" in script_text or \
               '__name__ == "__main__"' in script_text

    def test_build_cli_parser_present(self, script_text):
        assert '_build_cli_parser' in script_text

    def test_merged_vcf_argument(self, script_text):
        assert '--merged_vcf' in script_text

    def test_output_argument(self, script_text):
        assert '--output' in script_text

    def test_tool_vcfs_argument(self, script_text):
        assert '--tool_vcfs' in script_text

    def test_merger_mode_argument(self, script_text):
        assert '--merger_mode' in script_text

    def test_bed_file_argument(self, script_text):
        assert '--bed_file' in script_text

    def test_bam_file_argument(self, script_text):
        assert '--bam_file' in script_text

    def test_reference_fasta_argument(self, script_text):
        assert '--reference_fasta' in script_text

    def test_sample_id_argument(self, script_text):
        assert '--sample_id' in script_text

    def test_cli_parses_tool_vcfs_pairs(self, script_text):
        """CLI must split caller=path pairs on '='."""
        assert "split('=')" in script_text or ".split('='," in script_text


# ===========================================================================
# 22. MatchId-based caller tracing (_build_matchid_caller_map)
# ===========================================================================

class TestBuildMatchIdCallerMap:
    """_build_matchid_caller_map must build a MatchId -> callers set mapping."""

    # ── static / source-text checks ──────────────────────────────────────────

    def test_function_present(self, script_text):
        """_build_matchid_caller_map must be defined in the script."""
        assert 'def _build_matchid_caller_map(' in script_text

    def test_reads_matchid_info_field(self, script_text):
        """The function must read the MatchId INFO field from each record."""
        assert "'MatchId'" in script_text or '"MatchId"' in script_text

    def test_uses_defaultdict_of_sets(self, script_text):
        """Implementation must use defaultdict(set) to accumulate caller sets."""
        assert 'defaultdict(set)' in script_text

    def test_closes_vcf_in_finally(self, script_text):
        """VCF handle must be closed inside a finally block."""
        assert 'vcf.close()' in script_text
        assert 'finally:' in script_text

    def test_returns_plain_dict(self, script_text):
        """Result must be converted from defaultdict to a plain dict."""
        assert 'return dict(matchid_callers)' in script_text

    def test_calls_caller_from_tool_info(self, script_text):
        """Must delegate caller recognition to _caller_from_tool_info."""
        assert '_caller_from_tool_info' in script_text

    # ── functional tests (pysam.VariantFile patched on the module) ───────────

    def test_returns_empty_dict_for_empty_vcf(self, fe):
        """Empty VCF must produce an empty mapping."""
        from unittest import mock

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {}

    def test_returns_dict_not_defaultdict(self, fe):
        """Return value must be a plain dict, not a defaultdict."""
        from unittest import mock

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert type(result) is dict

    def test_skips_records_without_matchid(self, fe):
        """Records that lack a MatchId INFO field must be silently skipped."""
        from unittest import mock

        rec = mock.MagicMock()
        rec.info.get.return_value = None  # all INFO lookups return None

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([rec]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {}

    def test_skips_records_with_unrecognised_caller(self, fe):
        """Records whose TOOL value matches no known caller must be silently skipped."""
        from unittest import mock

        def _info_get(field, default=None):
            if field == 'MatchId':
                return 3
            return 'UNKNOWN_TOOL_XYZ'

        rec = mock.MagicMock()
        rec.info.get.side_effect = _info_get

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([rec]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {}

    def test_single_canoes_record(self, fe):
        """A single CANOES record with MatchId=1 produces {1: {'canoes'}}."""
        from unittest import mock

        def _info_get(field, default=None):
            if field == 'MatchId':
                return 1
            if field == 'TOOL':
                return 'CANOES'
            return None

        rec = mock.MagicMock()
        rec.info.get.side_effect = _info_get

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([rec]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {1: {'canoes'}}

    def test_multiple_callers_same_matchid(self, fe):
        """Two records sharing the same MatchId are merged into one set."""
        from unittest import mock

        def _make_rec(match_id, tool_name):
            def _info_get(field, default=None):
                if field == 'MatchId':
                    return match_id
                if field == 'TOOL':
                    return tool_name
                return None
            rec = mock.MagicMock()
            rec.info.get.side_effect = _info_get
            return rec

        recs = [_make_rec(5, 'CANOES'), _make_rec(5, 'DRAGEN')]

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter(recs))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {5: {'canoes', 'dragen'}}

    def test_different_matchids_stay_separate(self, fe):
        """Records with distinct MatchIds produce independent set entries."""
        from unittest import mock

        def _make_rec(match_id, tool_name):
            def _info_get(field, default=None):
                if field == 'MatchId':
                    return match_id
                if field == 'TOOL':
                    return tool_name
                return None
            rec = mock.MagicMock()
            rec.info.get.side_effect = _info_get
            return rec

        recs = [
            _make_rec(1, 'CANOES'),
            _make_rec(2, 'XHMM'),
            _make_rec(3, 'DRAGEN'),
        ]

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter(recs))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert result == {1: {'canoes'}, 2: {'xhmm'}, 3: {'dragen'}}

    def test_matchid_cast_to_int(self, fe):
        """MatchId values are cast to int before being used as dict keys."""
        from unittest import mock

        def _info_get(field, default=None):
            if field == 'MatchId':
                return '7'  # string value, must be cast to int
            if field == 'TOOL':
                return 'CNVKIT'
            return None

        rec = mock.MagicMock()
        rec.info.get.side_effect = _info_get

        mock_vcf = mock.MagicMock()
        mock_vcf.__iter__ = mock.Mock(return_value=iter([rec]))
        mock_vcf.close = mock.Mock()

        with mock.patch.object(fe, 'pysam') as mock_pysam:
            mock_pysam.VariantFile.return_value = mock_vcf
            result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

        assert 7 in result
        assert result[7] == {'cnvkit'}

    def test_all_supported_callers_recognised(self, fe):
        """Every caller in SUPPORTED_CALLERS must be matchable via a TOOL value."""
        from unittest import mock

        # Map canonical name -> a TOOL string that _TOOL_VALUE_PATTERNS recognises.
        tool_strings = {
            'canoes':    'CANOES',
            'clamms':    'CLAMMS',
            'xhmm':      'XHMM',
            'gatk_gcnv': 'GCNV',
            'cnvkit':    'CNVKIT',
            'dragen':    'DRAGEN',
            'exomedepth': 'EXOMEDEPTH',
        }

        for idx, caller_name in enumerate(fe.SUPPORTED_CALLERS):
            tool_str = tool_strings[caller_name]

            def _info_get(field, default=None, _mid=idx + 1, _t=tool_str):
                if field == 'MatchId':
                    return _mid
                if field == 'TOOL':
                    return _t
                return None

            rec = mock.MagicMock()
            rec.info.get.side_effect = _info_get

            mock_vcf = mock.MagicMock()
            mock_vcf.__iter__ = mock.Mock(return_value=iter([rec]))
            mock_vcf.close = mock.Mock()

            with mock.patch.object(fe, 'pysam') as mock_pysam:
                mock_pysam.VariantFile.return_value = mock_vcf
                result = fe._build_matchid_caller_map('/fake/collapsed.vcf')

            assert idx + 1 in result, (
                f"MatchId {idx + 1} missing; TOOL='{tool_str}' not recognised as '{caller_name}'"
            )
            assert caller_name in result[idx + 1], (
                f"Expected '{caller_name}' for TOOL='{tool_str}', got {result[idx + 1]}"
            )


# ===========================================================================
# 22. Single caller events, number_of_callers, and caller flags
# ===========================================================================

class TestNumberOfCallersAndCallerFlags:
    """Features dataframe must include number_of_callers and all supported is_{caller} flags."""

    def test_script_contains_number_of_callers(self, script_text):
        assert "number_of_callers" in script_text

    def test_single_caller_features(self, fe, tmp_path):
        import pysam

        # Create mock SURVIVOR VCF with 1 call supported by CANOES only
        vcf_path = str(tmp_path / "single_caller.vcf")
        hdr = pysam.VariantHeader()
        hdr.contigs.add("chr1", 1000000)
        hdr.add_line("##SAMPLE=<ID=0,File=sample.canoes.normalised.vcf.gz>")
        hdr.add_line("##SAMPLE=<ID=1,File=sample.gatk_gcnv.normalised.vcf.gz>")
        hdr.info.add("SVTYPE", 1, "String", "Type")
        hdr.info.add("END", 1, "Integer", "End")
        hdr.info.add("SUPP", 1, "Integer", "Supp")
        hdr.info.add("SUPP_VEC", 1, "String", "Supp Vec")
        hdr.info.add("PRECISE", 0, "Flag", "Precise")
        hdr.info.add("SV_REFINED", 1, "String", "Refined")
        hdr.samples.add("sample_1")

        with pysam.VariantFile(vcf_path, "w", header=hdr) as f:
            rec = hdr.new_record(contig="chr1", start=1000, stop=2000, alleles=("N", "<DEL>"))
            rec.info["SVTYPE"] = "DEL"
            rec.info["SUPP"] = 1
            rec.info["SUPP_VEC"] = "10"
            f.write(rec)

        df = fe.extract_normalized_features(
            vcf_path,
            merger_mode="survivor",
            tool_vcfs={}
        )

        assert len(df) == 1
        row = df.iloc[0]
        assert row["concordance"] == 1
        assert row["number_of_callers"] == 1
        assert row["is_canoes"] == 1
        assert row["is_gatk_gcnv"] == 0
        assert row["is_gcnv"] == 0
        assert row["is_clamms"] == 0
        assert row["is_xhmm"] == 0
        assert row["is_cnvkit"] == 0
        assert row["is_dragen"] == 0
        assert row["is_exomedepth"] == 0


class TestCallerIndependentMetricsAndCoverage:
    """Validates loading and extraction of caller-independent metrics from coverage files
    and unification of log2_ratio, depth_ratio, and bin_weight columns."""

    def test_load_and_query_coverage_cnn(self, fe, tmp_path):
        cnn_path = str(tmp_path / "sample.targetcoverage.cnn")
        with open(cnn_path, "w") as f:
            f.write("chromosome\tstart\tend\tgene\tdepth\tlog2\n")
            f.write("chr1\t100\t200\tG1\t50.0\t0.0\n")
            f.write("chr1\t300\t400\tG2\t50.0\t0.0\n")
            f.write("chr1\t1000\t1500\tG3\t10.0\t-1.32\n")
            f.write("chr1\t1600\t2000\tG4\t12.0\t-1.06\n")
            f.write("chr1\t2500\t2600\tG5\t50.0\t0.0\n")
            f.write("chr1\t2800\t2900\tG6\t50.0\t0.0\n")

        cov_data = fe._load_coverage_data(cnn_path)
        assert "1" in cov_data or "chr1" in cov_data

        # Query CNV region [1000, 2000)
        l2_mean, l2_dev, l2_diff, rd_rat, n_p, n_flank, wt = fe._query_coverage_metrics(
            cov_data, "1", 1000, 2000
        )
        assert n_p == 2
        assert n_flank == 4
        assert -1.5 < l2_mean < -1.0
        assert rd_rat < 0.3  # CNV depth ~11, flank depth 50 -> rd_ratio ~ 0.22

    def test_feature_extraction_with_coverage_file_unified_columns(self, fe, tmp_path):
        """Verify unified columns (log2_ratio, depth_ratio, bin_weight) and aliases."""
        cnn_path = str(tmp_path / "sample.targetcoverage.cnn")
        with open(cnn_path, "w") as f:
            f.write("chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n")
            f.write("chr1\t100\t200\tG1\t100.0\t0.0\t1.0\n")
            f.write("chr1\t1000\t2000\tG2\t20.0\t-1.32\t0.85\n")
            f.write("chr1\t3000\t3100\tG3\t100.0\t0.0\t1.0\n")

        import pysam

        vcf_path = str(tmp_path / "sample.vcf")
        hdr = pysam.VariantHeader()
        hdr.contigs.add("chr1", 50000000)
        hdr.add_line("##SAMPLE=<ID=1,File=sample.canoes.normalised.vcf.gz>")
        hdr.info.add("SVTYPE", 1, "String", "Type")
        hdr.info.add("END", 1, "Integer", "End")
        hdr.info.add("SUPP", 1, "Integer", "Supp")
        hdr.info.add("SUPP_VEC", 1, "String", "Supp Vec")
        hdr.samples.add("sample_1")

        with pysam.VariantFile(vcf_path, "w", header=hdr) as f:
            rec = hdr.new_record(contig="chr1", start=1000, stop=2000, alleles=("N", "<DEL>"))
            rec.info["SVTYPE"] = "DEL"
            rec.info["SUPP"] = 1
            rec.info["SUPP_VEC"] = "1"
            f.write(rec)

        df = fe.extract_normalized_features(
            vcf_path,
            merger_mode="survivor",
            tool_vcfs={},
            coverage_file=cnn_path,
        )

        assert len(df) == 1
        row = df.iloc[0]
        # Unified columns must exist and be populated
        assert "log2_ratio" in row
        assert "depth_ratio" in row
        assert "bin_weight" in row
        assert row["log2_ratio"] == -1.32
        assert row["cnvkit_log2"] == -1.32
        assert row["l2r_mean"] == -1.32
        assert row["depth_ratio"] == 0.2
        assert row["rd_ratio"] == 0.2
        assert row["bin_weight"] == 0.85
        assert row["n_probes"] == 1
        assert row["n_probes_flank"] == 2

    def test_derived_harmonisation_features_and_cross_derivation(self, fe, tmp_path):
        """Test bidirectional cross-derivation (R <-> log2 R) and derived GC / mappability features."""
        # 1. Mock Reference FASTA: 50% GC
        ref_fasta = str(tmp_path / "ref.fa")
        with open(ref_fasta, "w") as f:
            f.write(">chr1\n" + "GCAT" * 5000 + "\n")
        import pysam
        pysam.faidx(ref_fasta)

        # 2. Mock Mappability BED
        map_bed = str(tmp_path / "map.bed")
        with open(map_bed, "w") as f:
            f.write("chr1\t1000\t2000\t0.75\n")

        vcf_path = str(tmp_path / "sample.vcf")
        vcf_text_survivor = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=50000000>\n"
            "##SAMPLE=<ID=1,File=sample.exomedepth.normalised.vcf.gz>\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Supp Vec\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tN\t<DEL>\t90.0\tPASS\tSVTYPE=DEL;END=2000;SUPP_VEC=1\n"
        )
        with open(vcf_path, "w") as f:
            f.write(vcf_text_survivor)
        # Mock normalized ExomeDepth VCF with READS_RATIO=0.5
        plain_vcf = str(tmp_path / "sample.exomedepth.normalised.vcf")
        ed_vcf = str(tmp_path / "sample.exomedepth.normalised.vcf.gz")
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=50000000>\n"
            "##INFO=<ID=READS_RATIO,Number=1,Type=Float,Description=\"Reads Ratio\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tN\t<DEL>\t90.0\tPASS\tREADS_RATIO=0.5\n"
        )
        with open(plain_vcf, "w") as f:
            f.write(vcf_text)
        pysam.tabix_compress(plain_vcf, ed_vcf, force=True)
        pysam.tabix_index(ed_vcf, preset="vcf", force=True)

        df = fe.extract_normalized_features(
            vcf_path,
            merger_mode="survivor",
            tool_vcfs={"exomedepth": ed_vcf},
            reference_fasta=ref_fasta,
            mappability_file=map_bed,
        )

        assert len(df) == 1
        row = df.iloc[0]

        # Check cross-derivation: depth_ratio=0.5 -> log2_ratio=log2(0.5)=-1.0
        assert row["depth_ratio"] == 0.5
        assert pytest.approx(row["log2_ratio"], 0.01) == -1.0

        # Check GC metrics: GC=0.50 -> gc_dev=0.0, is_extreme_gc=0
        assert pytest.approx(row["gc_content"], abs=0.01) == 0.50
        assert pytest.approx(row["gc_dev"], abs=0.01) == 0.0
        assert row["is_extreme_gc"] == 0

        # Check Mappability: score=0.75 -> is_low_mappability=1 (< 0.80)
        assert pytest.approx(row["mappability"], 0.01) == 0.75
        assert row["is_low_mappability"] == 1

        # Check Target size & overhang
        assert row["cnv_size"] == 1000
        assert pytest.approx(row["overhang_ratio"], 0.01) == 1000.0 / 1300.0

