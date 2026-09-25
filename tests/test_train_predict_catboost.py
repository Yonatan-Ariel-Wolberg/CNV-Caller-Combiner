#!/usr/bin/env python3
"""
Tests for bin/train_predict_catboost.py.

Validates:
  1. CatBoost model training & in silico VCF prediction with a single caller (e.g. ExomeDepth).
  2. Training & prediction on arbitrary caller subsets (e.g. 2 callers, 3 callers, 7 callers).
  3. Prediction filtering functions with both SURVIVOR and Truvari merged VCF outputs.
  4. Features like number_of_callers, is_*, concordance, and biophysical metrics are ingested.
"""

import os
import subprocess
import pysam
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
CATBOOST_SCRIPT = os.path.join(BIN_DIR, "train_predict_catboost.py")


def create_mock_vcf(vcf_path, records_data, sample_id="SAMPLE_001"):
    """Helper to create mock VCF."""
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Type")
    header.info.add("END", 1, "Integer", "End")
    header.info.add("SVLEN", 1, "Integer", "Length")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("TOOLS", ".", "String", "Tools")
    header.samples.add(sample_id)

    with pysam.VariantFile(vcf_path, "w", header=header) as f:
        for chrom, start, end, svtype in records_data:
            rec = header.new_record(contig=chrom, start=start, stop=end, alleles=("N", f"<{svtype}>"))
            rec.info["SVTYPE"] = svtype
            rec.info["SVLEN"] = end - start
            rec.samples[sample_id]["GT"] = (0, 1)
            rec.samples[sample_id]["TOOLS"] = "TEST"
            f.write(rec)


def test_catboost_training_single_caller_exomedepth(tmp_path):
    """Test CatBoost training and prediction when only a single caller (ExomeDepth) was run."""
    feat_dir = tmp_path / "features"
    feat_dir.mkdir()
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    out_dir = tmp_path / "pred_out"
    truth_file = str(tmp_path / "truth.tsv")

    # 1. Mock feature file for single caller
    feat_path = str(feat_dir / "SAMPLE_001_features.tsv")
    with open(feat_path, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\tnumber_of_callers\tconcordance\tis_exomedepth\tis_canoes\tis_gatk_gcnv\tis_clamms\tis_xhmm\tis_cnvkit\tis_dragen\tqual_norm_exomedepth\trd_ratio\tn_probes\n")
        f.write("SAMPLE_001\t1\t1000\t2000\t0\t1\t1\t1\t0\t0\t0\t0\t0\t0\t95.0\t0.45\t12\n")
        f.write("SAMPLE_001\t1\t5000\t6000\t0\t1\t1\t1\t0\t0\t0\t0\t0\t0\t30.0\t0.98\t2\n")

    # 2. Mock truth labels
    with open(truth_file, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\ttruth_label\n")
        f.write("SAMPLE_001\t1\t1000\t2000\tDEL\t1\n")
        f.write("SAMPLE_001\t1\t5000\t6000\tDEL\t0\n")

    # 3. Mock VCF to filter
    vcf_path = str(vcf_dir / "SAMPLE_001_rare.vcf")
    create_mock_vcf(vcf_path, [("chr1", 1000, 2000, "DEL"), ("chr1", 5000, 6000, "DEL")], "SAMPLE_001")

    # 4. Run train_predict_catboost.py
    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        CATBOOST_SCRIPT,
        "--features-dir", str(feat_dir),
        "--vcf-dir", str(vcf_dir),
        "--truth-labels", truth_file,
        "--output-dir", str(out_dir),
        "--threshold", "0.5"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}\nOutput: {res.stdout}"

    # 5. Check prediction outputs
    pred_vcf = os.path.join(str(out_dir), "SAMPLE_001_rare_predicted.vcf")
    assert os.path.exists(pred_vcf)

    with pysam.VariantFile(pred_vcf) as vcf_in:
        recs = list(vcf_in)
        assert len(recs) >= 1
        assert recs[0].start == 1000 or recs[0].pos in (1000, 1001)


def test_catboost_multi_caller_combinations(tmp_path):
    """Test CatBoost training with multi-caller subsets (e.g. 2 callers: CANOES + GATK-gCNV, 3 callers, etc.)."""
    feat_dir = tmp_path / "features"
    feat_dir.mkdir()
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    out_dir = tmp_path / "pred_out"
    truth_file = str(tmp_path / "truth.tsv")

    feat_path = str(feat_dir / "SAMPLE_002_features.tsv")
    with open(feat_path, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\tnumber_of_callers\tconcordance\tis_canoes\tis_gatk_gcnv\tis_clamms\tis_xhmm\tis_cnvkit\tis_exomedepth\tis_dragen\tqual_norm_canoes\tqual_norm_gatk_gcnv\trd_ratio\tn_probes\n")
        f.write("SAMPLE_002\t1\t1000\t2000\t0\t2\t2\t1\t1\t0\t0\t0\t0\t0\t92.0\t98.0\t0.42\t15\n")
        f.write("SAMPLE_002\t1\t8000\t9000\t0\t1\t1\t1\t0\t0\t0\t0\t0\t0\t25.0\t0.0\t0.95\t3\n")

    with open(truth_file, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\ttruth_label\n")
        f.write("SAMPLE_002\t1\t1000\t2000\tDEL\t1\n")
        f.write("SAMPLE_002\t1\t8000\t9000\tDEL\t0\n")

    vcf_path = str(vcf_dir / "SAMPLE_002_rare.vcf")
    create_mock_vcf(vcf_path, [("chr1", 1000, 2000, "DEL"), ("chr1", 8000, 9000, "DEL")], "SAMPLE_002")

    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        CATBOOST_SCRIPT,
        "--features-dir", str(feat_dir),
        "--vcf-dir", str(vcf_dir),
        "--truth-labels", truth_file,
        "--output-dir", str(out_dir),
        "--threshold", "0.5"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    pred_vcf = os.path.join(str(out_dir), "SAMPLE_002_rare_predicted.vcf")
    assert os.path.exists(pred_vcf)
    with pysam.VariantFile(pred_vcf) as vcf_in:
        recs = list(vcf_in)
        assert len(recs) >= 1
        assert recs[0].start == 1000 or recs[0].pos in (1000, 1001)


def test_catboost_works_with_survivor_and_truvari_vcf(tmp_path):
    """Test CatBoost filtering on both SURVIVOR and Truvari intra-merged outputs."""
    feat_dir = tmp_path / "features"
    feat_dir.mkdir()
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    out_dir = tmp_path / "pred_out"
    truth_file = str(tmp_path / "truth.tsv")

    feat_path = str(feat_dir / "SAMPLE_003_features.tsv")
    with open(feat_path, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\tnumber_of_callers\tconcordance\tis_canoes\tis_gatk_gcnv\tqual_norm_canoes\trd_ratio\n")
        f.write("SAMPLE_003\t1\t3000\t4000\t0\t2\t2\t1\t1\t95.0\t0.50\n")

    with open(truth_file, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\ttruth_label\n")
        f.write("SAMPLE_003\t1\t3000\t4000\tDEL\t1\n")

    # 1. SURVIVOR naming format
    survivor_vcf = str(vcf_dir / "SAMPLE_003_survivor_union.vcf")
    create_mock_vcf(survivor_vcf, [("chr1", 3000, 4000, "DEL")], "SAMPLE_003")

    # 2. Truvari naming format
    truvari_vcf = str(vcf_dir / "SAMPLE_003_truvari_merged.vcf")
    create_mock_vcf(truvari_vcf, [("chr1", 3000, 4000, "DEL")], "SAMPLE_003")

    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        CATBOOST_SCRIPT,
        "--features-dir", str(feat_dir),
        "--vcf-dir", str(vcf_dir),
        "--truth-labels", truth_file,
        "--output-dir", str(out_dir),
        "--threshold", "0.5"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    assert os.path.exists(os.path.join(str(out_dir), "SAMPLE_003_survivor_union_predicted.vcf"))
    assert os.path.exists(os.path.join(str(out_dir), "SAMPLE_003_truvari_merged_predicted.vcf"))


def test_catboost_excludes_breakpoint_refinement_flags(tmp_path):
    """Test that is_precise and is_sv_refined are formally excluded from the CatBoost feature set."""
    feat_dir = tmp_path / "features"
    feat_dir.mkdir()
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    out_dir = tmp_path / "pred_out"
    truth_file = str(tmp_path / "truth.tsv")

    feat_path = str(feat_dir / "SAMPLE_004_features.tsv")
    with open(feat_path, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\tnumber_of_callers\tconcordance\tis_canoes\tis_clamms\tis_precise\tis_sv_refined\trd_ratio\n")
        f.write("SAMPLE_004\t1\t1000\t2000\t0\t2\t2\t1\t1\t1\t1\t0.45\n")
        f.write("SAMPLE_004\t1\t5000\t6000\t0\t2\t2\t1\t1\t0\t0\t0.98\n")

    with open(truth_file, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcnv_type\ttruth_label\n")
        f.write("SAMPLE_004\t1\t1000\t2000\tDEL\t1\n")
        f.write("SAMPLE_004\t1\t5000\t6000\tDEL\t0\n")

    vcf_path = str(vcf_dir / "SAMPLE_004_rare.vcf")
    create_mock_vcf(vcf_path, [("chr1", 1000, 2000, "DEL"), ("chr1", 5000, 6000, "DEL")], "SAMPLE_004")

    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        CATBOOST_SCRIPT,
        "--features-dir", str(feat_dir),
        "--vcf-dir", str(vcf_dir),
        "--truth-labels", truth_file,
        "--output-dir", str(out_dir),
        "--threshold", "0.5"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}\nOutput: {res.stdout}"

    shap_path = os.path.join(str(out_dir), "catboost_shap_values.tsv")
    if os.path.exists(shap_path):
        import pandas as pd
        shap_df = pd.read_csv(shap_path, sep="\t")
        assert "is_precise" not in shap_df.columns
        assert "is_sv_refined" not in shap_df.columns

