#!/usr/bin/env python3
"""
Unit tests for the statistical significance and comparative benchmark engine
(bin/compare_pipeline_benchmarks.py).
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import pytest

# Add bin/ to import path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'bin'))

from compare_pipeline_benchmarks import (
    compute_f_beta,
    compute_mcc,
    compute_odds_ratio,
    compute_relative_risk_sensitivity,
    mcnemar_test,
    cochrans_q_test,
    benjamini_hochberg,
    pages_l_test,
    calculate_nri_idi,
    cluster_bootstrap_metric_deltas,
    PipelineBenchmarkComparator
)


def test_compute_f_beta():
    # Perfect precision and recall -> F2 should be 1.0
    assert compute_f_beta(1.0, 1.0, 2.0) == 1.0
    # Zero precision or recall -> 0.0
    assert compute_f_beta(0.0, 1.0, 2.0) == 0.0
    assert compute_f_beta(1.0, 0.0, 2.0) == 0.0
    # F2 puts more weight on recall than precision:
    # High recall (0.9), low precision (0.5) -> F2 should be higher than F1
    f2_val = compute_f_beta(0.5, 0.9, 2.0)
    assert 0.7 < f2_val < 0.9


def test_compute_mcc():
    # Perfect classification
    assert compute_mcc(100, 0, 0, 100) == 1.0
    # Complete inversion
    assert compute_mcc(0, 100, 100, 0) == -1.0
    # Random / orthogonal
    assert compute_mcc(50, 50, 50, 50) == 0.0


def test_compute_odds_ratio():
    or_val, ci_low, ci_high = compute_odds_ratio(80, 20, 10, 90, 60, 40, 30, 70)
    assert or_val > 1.0
    assert ci_low < or_val < ci_high


def test_compute_relative_risk_sensitivity():
    rr, ci_low, ci_high = compute_relative_risk_sensitivity(90, 10, 60, 40)
    assert rr > 1.0
    assert ci_low < rr < ci_high


def test_mcnemar_test():
    # Equal discordance -> chi2 near 0, p near 1
    chi2, p = mcnemar_test(10, 10)
    assert p == 1.0
    # Asymmetric discordance -> significant p
    chi2, p = mcnemar_test(30, 2)
    assert p < 0.001


def test_cochrans_q_test():
    # Matrix of 20 samples across 4 conditions
    rng = np.random.default_rng(42)
    mat = np.column_stack([
        rng.binomial(1, 0.9, 30),
        rng.binomial(1, 0.7, 30),
        rng.binomial(1, 0.5, 30),
        rng.binomial(1, 0.3, 30)
    ])
    q_stat, p = cochrans_q_test(mat)
    assert q_stat > 0
    assert p < 0.05


def test_benjamini_hochberg():
    p_vals = [0.001, 0.01, 0.04, 0.20, 0.80]
    adj = benjamini_hochberg(p_vals)
    assert len(adj) == len(p_vals)
    # Monotonicity check
    assert adj[0] <= adj[1] <= adj[2] <= adj[3] <= adj[4]


def test_pages_l_test():
    # Strictly increasing ranks
    matrix = np.array([
        [1.0, 2.0, 3.0, 4.0],
        [1.1, 2.2, 3.1, 4.2],
        [0.9, 2.1, 2.9, 3.8],
        [1.2, 2.0, 3.2, 4.1]
    ])
    l_stat, p = pages_l_test(matrix)
    assert l_stat > 0
    assert p < 0.05


def test_calculate_nri_idi():
    rng = np.random.default_rng(42)
    y_true = rng.binomial(1, 0.5, 100)
    p_base = rng.uniform(0.2, 0.6, 100)
    # Better model gives higher probs for events
    p_new = np.where(y_true == 1, p_base + 0.2, p_base - 0.2)
    p_new = np.clip(p_new, 0.0, 1.0)
    nri, p_nri, idi, p_idi = calculate_nri_idi(y_true, p_base, p_new)
    assert nri > 0
    assert idi > 0


def test_cluster_bootstrap_metric_deltas():
    rng = np.random.default_rng(42)
    samples_a = {
        f"sample_{i}": pd.DataFrame({'TP': [rng.integers(15, 25)], 'FP': [rng.integers(1, 5)], 'FN': [rng.integers(1, 5)], 'TN': [rng.integers(80, 100)]})
        for i in range(10)
    }
    samples_b = {
        f"sample_{i}": pd.DataFrame({'TP': [rng.integers(10, 15)], 'FP': [rng.integers(5, 10)], 'FN': [rng.integers(5, 10)], 'TN': [rng.integers(70, 90)]})
        for i in range(10)
    }
    deltas = cluster_bootstrap_metric_deltas(samples_a, samples_b, n_resamples=50, seed=42)
    assert 'delta_sensitivity' in deltas
    assert 'delta_precision' in deltas
    assert 'delta_f2' in deltas
    assert 'delta_mcc' in deltas
    # Sensitivity A > B so delta should be positive
    val, low, high = deltas['delta_sensitivity']
    assert val > 0


def test_pipeline_benchmark_comparator_end_to_end(tmp_path):
    comparator = PipelineBenchmarkComparator(seed=42)
    df = comparator.run_all_dimensions()
    assert len(df) > 0
    assert 'comparison_dimension' in df.columns
    assert 'p_adjusted_fdr' in df.columns
    assert 'odds_ratio' in df.columns

    # Test CLI execution
    tsv_out = tmp_path / "comparisons.tsv"
    report_out = tmp_path / "report.md"
    cmd = [
        sys.executable,
        os.path.join(REPO_ROOT, 'bin', 'compare_pipeline_benchmarks.py'),
        '--output-tsv', str(tsv_out),
        '--output-report', str(report_out),
        '--bootstrap-iterations', '50'
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert tsv_out.exists()
    assert report_out.exists()
