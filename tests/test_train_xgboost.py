#!/usr/bin/env python3
"""
Tests for bin/train_xgboost.py.

Validates:
  1. SUPPORTED_CALLERS lists all seven pipeline callers.
  2. MIN_CALLERS_FOR_TRAINING equals 2 (concordance requires >= 2 callers).
  3. validate_min_callers raises ValueError for 0 or 1 caller and passes for >= 2.
  4. prepare_training_data uses named is_* columns when present (feature-
     extraction output) and validates minimum callers.
  5. prepare_training_data falls back to the legacy supp_vec path and
     validates minimum callers there too.
"""

import importlib.util
import os
import re
import sys
import ast

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Locate and import the script
# ---------------------------------------------------------------------------

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT_PATH = os.path.join(BIN_DIR, 'train_xgboost.py')


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _import_module():
    spec = importlib.util.spec_from_file_location('train_xgboost', SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def tx():
    """Return the train_xgboost module."""
    return _import_module()

EXPECTED_TRUTH_LABEL_COLUMNS = ['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label']
# merge_features_with_truth_labels joins on cnv_type_norm where cnv_type values
# are normalized into canonical labels (e.g. 1/DUP -> DUP, 0/DEL -> DEL).
EXPECTED_MERGE_KEYS = ['sample_id', 'chrom', 'start', 'end', 'cnv_type_norm']


# ===========================================================================
# 1. SUPPORTED_CALLERS constant
# ===========================================================================

class TestSupportedCallers:

    def test_supported_callers_present(self, script_text):
        assert 'SUPPORTED_CALLERS' in script_text

    def test_supported_callers_has_seven_entries(self, tx):
        assert len(tx.SUPPORTED_CALLERS) == 7

    def test_supported_callers_contains_canoes(self, tx):
        assert 'canoes' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_clamms(self, tx):
        assert 'clamms' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_xhmm(self, tx):
        assert 'xhmm' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_gatk_gcnv(self, tx):
        assert 'gatk_gcnv' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_cnvkit(self, tx):
        assert 'cnvkit' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_dragen(self, tx):
        assert 'dragen' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_indelible(self, tx):
        assert 'indelible' in tx.SUPPORTED_CALLERS


# ===========================================================================
# 2. MIN_CALLERS_FOR_TRAINING constant
# ===========================================================================

class TestMinCallersConstant:

    def test_min_callers_present(self, script_text):
        assert 'MIN_CALLERS_FOR_TRAINING' in script_text

    def test_min_callers_equals_two(self, tx):
        assert tx.MIN_CALLERS_FOR_TRAINING == 2

    def test_min_callers_less_than_supported_callers(self, tx):
        assert tx.MIN_CALLERS_FOR_TRAINING < len(tx.SUPPORTED_CALLERS)


# ===========================================================================
# 3. validate_min_callers function
# ===========================================================================

class TestValidateMinCallers:

    def test_validate_min_callers_present(self, script_text):
        assert 'def validate_min_callers(' in script_text

    def test_zero_callers_raises(self, tx):
        with pytest.raises(ValueError, match='At least 2'):
            tx.validate_min_callers([])

    def test_one_caller_raises(self, tx):
        with pytest.raises(ValueError, match='At least 2'):
            tx.validate_min_callers(['is_canoes'])

    def test_two_callers_passes(self, tx):
        tx.validate_min_callers(['is_canoes', 'is_clamms'])  # must not raise

    def test_seven_callers_passes(self, tx):
        cols = [f'is_{c}' for c in tx.SUPPORTED_CALLERS]
        tx.validate_min_callers(cols)  # must not raise

    def test_error_message_mentions_supported_callers(self, tx):
        with pytest.raises(ValueError, match='Supported callers'):
            tx.validate_min_callers([])

    def test_error_message_mentions_found_count(self, tx):
        with pytest.raises(ValueError, match='Found 1'):
            tx.validate_min_callers(['is_canoes'])


# ===========================================================================
# 4. prepare_training_data – named is_* columns (feature-extraction output)
# ===========================================================================

class TestPrepareTrainingDataNamedColumns:

    def _make_df_named(self, callers=None):
        """Build a minimal feature-extraction-style DataFrame."""
        if callers is None:
            callers = ['canoes', 'clamms']
        data = {
            'concordance': [1, 2],
            'cnv_size': [5000, 10000],
        }
        for c in callers:
            data[f'is_{c}'] = [1, 0]
            data[f'qual_norm_{c}'] = [80.0, np.nan]
        return pd.DataFrame(data)

    def test_returns_dataframe(self, tx):
        df = self._make_df_named()
        result = tx.prepare_training_data(df)
        assert isinstance(result, pd.DataFrame)

    def test_is_columns_preserved(self, tx):
        df = self._make_df_named(['canoes', 'clamms', 'xhmm'])
        result = tx.prepare_training_data(df)
        assert 'is_canoes' in result.columns
        assert 'is_clamms' in result.columns
        assert 'is_xhmm' in result.columns

    def test_one_caller_raises(self, tx):
        df = self._make_df_named(['canoes'])
        with pytest.raises(ValueError):
            tx.prepare_training_data(df)

    def test_no_supp_vec_column_needed(self, tx):
        """Feature-extraction output has no supp_vec; must not raise KeyError."""
        df = self._make_df_named(['canoes', 'clamms'])
        assert 'supp_vec' not in df.columns
        result = tx.prepare_training_data(df)  # must not raise
        assert result is not None


# ===========================================================================
# 5. prepare_training_data – legacy supp_vec path
# ===========================================================================

class TestPrepareTrainingDataLegacy:

    def _make_df_supp(self, n_callers=7):
        """Build a minimal DataFrame with a supp_vec column."""
        vec = '1' * n_callers
        return pd.DataFrame({
            'supp_vec': [vec, vec],
            'concordance': [n_callers, n_callers],
        })

    def test_legacy_path_seven_callers(self, tx):
        df = self._make_df_supp(7)
        result = tx.prepare_training_data(df, num_callers=7)
        assert 'supp_vec' not in result.columns
        for i in range(7):
            assert f'caller_{i}_flag' in result.columns

    def test_legacy_path_two_callers(self, tx):
        df = self._make_df_supp(2)
        result = tx.prepare_training_data(df, num_callers=2)
        assert 'caller_0_flag' in result.columns
        assert 'caller_1_flag' in result.columns

    def test_legacy_path_one_caller_raises(self, tx):
        df = self._make_df_supp(1)
        with pytest.raises(ValueError):
            tx.prepare_training_data(df, num_callers=1)

    def test_legacy_path_zero_callers_raises(self, tx):
        df = self._make_df_supp(0)
        with pytest.raises(ValueError):
            tx.prepare_training_data(df, num_callers=0)


# ===========================================================================
# 6. truth-label schema / join contract
# ===========================================================================

class TestTruthLabelContract:

    def test_truth_labels_required_columns_include_cnv_type(self, script_text):
        module = ast.parse(script_text)
        required_cols = None
        main_fn = next(
            (n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'main'),
            None
        )
        assert main_fn is not None
        for node in ast.walk(main_fn):
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == 'required_cols':
                    required_cols = ast.literal_eval(node.value)
                    break
            if required_cols is not None:
                break
        assert isinstance(required_cols, set)
        for col in EXPECTED_TRUTH_LABEL_COLUMNS:
            assert col in required_cols

    def test_merge_uses_normalized_cnv_type_in_join_keys(self, script_text):
        module = ast.parse(script_text)
        join_keys = None
        merge_fn = next(
            (n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'merge_features_with_truth_labels'),
            None
        )
        assert merge_fn is not None
        for node in ast.walk(merge_fn):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                if node.func.attr == 'merge':
                    for kw in node.keywords:
                        if kw.arg == 'on':
                            join_keys = ast.literal_eval(kw.value)
                            break
            if join_keys is not None:
                break
        assert join_keys is not None
        assert join_keys == EXPECTED_MERGE_KEYS

    def test_truth_labels_cli_help_mentions_cnv_type(self, script_text):
        assert 'sample_id, chrom, start, end, cnv_type, truth_label' in script_text

    def test_cli_supports_probe_overlap_matching(self, script_text):
        assert '--probes_bed' in script_text
        assert '--min_shared_probes' in script_text

    def test_reports_probe_overlap_matches(self, script_text):
        assert 'Probe-overlap matches' in script_text

    def test_cli_supports_roc_and_pr_outputs(self, script_text):
        assert '--output_roc_plot' in script_text
        assert '--output_pr_plot' in script_text
        assert '--output_roc_data' in script_text
        assert '--output_pr_data' in script_text

    def test_cli_supports_shap_outputs(self, script_text):
        assert '--output_shap_values' in script_text
        assert '--output_shap_summary_plot' in script_text
        assert '--output_shap_beeswarm_plot' in script_text

    def test_reports_auc_metrics(self, script_text):
        assert 'ROC AUC:' in script_text
        assert 'PR AUC:' in script_text


# ===========================================================================
# 7. probe-overlap fallback matching
# ===========================================================================

class TestProbeOverlapMatching:

    def _features_df(self):
        return pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'chrom': ['chr1', 'chr1'],
            'start': [100, 500],
            'end': [200, 600],
            'cnv_type': [1, 0],
            'concordance': [2, 2],
            'is_canoes': [1, 1],
            'is_clamms': [1, 1],
        })

    def _labels_df(self):
        # S1 exact coordinate mismatch but same probes should match by fallback.
        # S2 has identical chr1:500-600 coordinates in both DataFrames and is
        # therefore counted as an exact-key match.
        return pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'chrom': ['chr1', 'chr1'],
            'start': [110, 500],
            'end': [190, 600],
            'cnv_type': ['DUP', 'DEL'],
            'truth_label': [1, 0],
        })

    def _write_probe_bed(self, tmp_path):
        bed = tmp_path / 'probes.bed'
        bed.write_text(
            "chr1\t50\t120\n"
            "chr1\t130\t170\n"
            "chr1\t180\t210\n"
            "chr1\t500\t550\n"
            "chr1\t550\t620\n"
        )
        return str(bed)

    def test_probe_overlap_fallback_labels_unmatched_exact_rows(self, tx, tmp_path):
        merged, exact_count, probe_count = tx.merge_features_with_truth_labels(
            self._features_df(),
            self._labels_df(),
            probes_bed=self._write_probe_bed(tmp_path),
            min_shared_probes=1,
        )
        assert len(merged) == 2
        assert exact_count == 1
        assert probe_count == 1

    def test_min_shared_probes_threshold_is_enforced(self, tx, tmp_path):
        merged, exact_count, probe_count = tx.merge_features_with_truth_labels(
            self._features_df(),
            self._labels_df(),
            probes_bed=self._write_probe_bed(tmp_path),
            min_shared_probes=4,
        )
        # only the exact S2 match remains when threshold is too strict for S1
        assert len(merged) == 1
        assert exact_count == 1
        assert probe_count == 0
