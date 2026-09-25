#!/usr/bin/env python3
"""
Unit tests for internal dataset allele frequency filtering (bin/filter-rare-cnvs-ddduk-90.py).

Verifies that:
1. Allele frequency filtering focuses strictly and internally on our dataset/cohort samples.
2. External databases (e.g. ClinVar, gnomAD, AnnotSV) are NOT used or consulted for AF calculation.
3. Overlapping CNVs across cohort samples are properly clustered (reciprocal overlap >= 50%).
4. Common variants exceeding max_af are filtered out, while rare variants and singletons are preserved.
5. Empty inputs and region exports behave robustly without breaking Nextflow channels.
"""

import os
import subprocess
import pandas as pd
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
FILTER_SCRIPT = os.path.join(BIN_DIR, "filter-rare-cnvs-ddduk-90.py")


def test_internal_allele_frequency_filtering(tmp_path):
    """Test that AF is computed strictly from internal dataset samples and filtered correctly."""
    # Cohort of 10 samples
    total_samples = 10
    max_af = 0.20  # Keep variants with count <= 2

    # Create mock cohort TSV (flattened from bcftools: CHROM, POS, END, SVTYPE, SAMPLE, GT)
    tsv_content = (
        # Variant 1: Common in cohort (3 samples: S1, S2, S3) -> AF = 3/10 = 0.30 -> Exceeds max_af
        "chr1\t10000\t20000\tDEL\tS1\t0/1\n"
        "chr1\t10050\t20050\tDEL\tS2\t0/1\n"
        "chr1\t10000\t20000\tDEL\tS3\t0/1\n"
        # Variant 2: Rare in cohort (1 sample: S4) -> AF = 1/10 = 0.10 -> Retained
        "chr1\t50000\t60000\tDUP\tS4\t0/1\n"
        # Variant 3: Rare in cohort (2 samples: S5, S6) -> AF = 2/10 = 0.20 -> Retained
        "chr2\t100000\t150000\tDEL\tS5\t0/1\n"
        "chr2\t100100\t150000\tDEL\tS6\t0/1\n"
    )

    in_tsv = tmp_path / "cohort_input.tsv"
    in_tsv.write_text(tsv_content)

    out_csv = tmp_path / "rare_cnvs.csv"
    out_regions = tmp_path / "rare_regions.tsv"

    cmd = [
        "python3", FILTER_SCRIPT,
        "--input-tsv", str(in_tsv),
        "--output-csv", str(out_csv),
        "--total-samples", str(total_samples),
        "--max-af", str(max_af),
        "--output-regions", str(out_regions)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Script failed: {result.stderr}"

    # Verify CSV results
    df = pd.read_csv(out_csv)
    assert len(df) == 3  # S4 (DUP), S5 (DEL), S6 (DEL)
    assert "S1" not in df["SAMPLE"].values
    assert "S2" not in df["SAMPLE"].values
    assert "S3" not in df["SAMPLE"].values

    # Check internal AF calculation
    s4_row = df[df["SAMPLE"] == "S4"].iloc[0]
    assert s4_row["COUNT"] == 1
    assert pytest.approx(s4_row["INTERNAL_AF"], 0.01) == 0.10

    s5_row = df[df["SAMPLE"] == "S5"].iloc[0]
    assert s5_row["COUNT"] == 2
    assert pytest.approx(s5_row["INTERNAL_AF"], 0.01) == 0.20

    # Verify regions file
    assert out_regions.exists()
    regions_lines = [line.strip() for line in out_regions.read_text().splitlines() if line.strip()]
    # Should contain unique rare regions
    assert len(regions_lines) == 3


def test_singletons_preserved_in_small_cohorts(tmp_path):
    """In small cohorts, singletons (count=1) should be preserved even if count/total > max_af."""
    total_samples = 2
    max_af = 0.05  # 2 * 0.05 = 0.1, but singleton min count rule sets max_allowed_count = 1

    tsv_content = (
        "chr1\t1000\t5000\tDEL\tS1\t0/1\n"
        "chr1\t1000\t5000\tDEL\tS2\t0/1\n"  # 2 samples -> common in cohort -> filtered
        "chr2\t2000\t6000\tDUP\tS1\t0/1\n"  # 1 sample -> singleton -> preserved
    )
    in_tsv = tmp_path / "cohort_small.tsv"
    in_tsv.write_text(tsv_content)
    out_csv = tmp_path / "rare_small.csv"

    cmd = [
        "python3", FILTER_SCRIPT,
        "--input-tsv", str(in_tsv),
        "--output-csv", str(out_csv),
        "--total-samples", str(total_samples),
        "--max-af", str(max_af)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0

    df = pd.read_csv(out_csv)
    assert len(df) == 1
    assert df.iloc[0]["SAMPLE"] == "S1"
    assert df.iloc[0]["CHROM"] == "chr2"
    assert df.iloc[0]["COUNT"] == 1


def test_external_databases_not_consulted(tmp_path):
    """Providing an external database like AnnotSV/ClinVar is ignored and does not affect internal AF."""
    total_samples = 5
    max_af = 0.20

    tsv_content = (
        "chr1\t1000\t5000\tDEL\tS1\t0/1\n"
    )
    in_tsv = tmp_path / "cohort_ext.tsv"
    in_tsv.write_text(tsv_content)

    dummy_annotsv = tmp_path / "dummy_annotsv.tsv"
    dummy_annotsv.write_text(
        "Annotation_mode\tSV_chrom\tSV_start\tSV_end\tSV_type\tACMG_class\n"
        "full\tchr1\t1000\t5000\tDEL\t5\n"
    )

    out_csv_no_ext = tmp_path / "no_ext.csv"
    out_csv_with_ext = tmp_path / "with_ext.csv"

    # Run without external DB
    subprocess.run([
        "python3", FILTER_SCRIPT,
        "--input-tsv", str(in_tsv),
        "--output-csv", str(out_csv_no_ext),
        "--total-samples", str(total_samples),
        "--max-af", str(max_af)
    ], check=True)

    # Run with external DB flag
    res = subprocess.run([
        "python3", FILTER_SCRIPT,
        "--input-tsv", str(in_tsv),
        "--output-csv", str(out_csv_with_ext),
        "--total-samples", str(total_samples),
        "--max-af", str(max_af),
        "--annotsv-tsv", str(dummy_annotsv)
    ], capture_output=True, text=True)

    assert res.returncode == 0
    assert "External annotations (--annotsv-tsv) are ignored during allele frequency filtering" in res.stderr

    df_no_ext = pd.read_csv(out_csv_no_ext)
    df_with_ext = pd.read_csv(out_csv_with_ext)

    # Output columns and content must be identical
    pd.testing.assert_frame_equal(df_no_ext, df_with_ext)
    assert "ACMG_class" not in df_with_ext.columns


def test_empty_input_handling(tmp_path):
    """Empty input TSV generates valid empty output CSV and regions file without failing."""
    empty_tsv = tmp_path / "empty.tsv"
    empty_tsv.write_text("")
    out_csv = tmp_path / "empty_out.csv"
    out_regions = tmp_path / "empty_regions.tsv"

    cmd = [
        "python3", FILTER_SCRIPT,
        "--input-tsv", str(empty_tsv),
        "--output-csv", str(out_csv),
        "--total-samples", "10",
        "--max-af", "0.05",
        "--output-regions", str(out_regions)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0

    df = pd.read_csv(out_csv)
    assert df.empty
    assert list(df.columns) == ['CHROM', 'POS', 'END', 'SVTYPE', 'COUNT', 'INTERNAL_AF', 'SAMPLE', 'GT']
    assert out_regions.exists()
    assert out_regions.read_text() == ""
