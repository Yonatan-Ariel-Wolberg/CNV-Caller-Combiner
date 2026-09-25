#!/usr/bin/env python3
"""
Tests for ExomeDepth integration into the pipeline, including:
  1. modules/callers/modules-exomedepth.nf structure, processes, and workflow.
  2. bin/exomedepth_csv_to_vcf.py conversion logic, headers, and fields.
  3. ML model and feature extraction integration.
"""

import os
import re
import sys
import tempfile
import importlib.util
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
NF_MODULE = os.path.join(REPO_ROOT, 'modules', 'callers', 'modules-exomedepth.nf')
CONVERTER = os.path.join(REPO_ROOT, 'bin', 'exomedepth_csv_to_vcf.py')


def _read_file(path):
    with open(path) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def nf_text():
    return _read_file(NF_MODULE)


@pytest.fixture(scope='module')
def converter_text():
    return _read_file(CONVERTER)


# ===========================================================================
# 1. Nextflow Module DSL2 & Structure
# ===========================================================================

class TestExomeDepthModuleStructure:

    def test_file_exists(self):
        assert os.path.isfile(NF_MODULE), "modules/callers/modules-exomedepth.nf must exist"

    def test_dsl2_enabled(self, nf_text):
        assert "nextflow.enable.dsl=2" in nf_text

    def test_workflow_defined(self, nf_text):
        assert "workflow EXOMEDEPTH {" in nf_text

    def test_workflow_emits(self, nf_text):
        assert "sorted_vcf" in nf_text
        assert "sorted_vcf_index" in nf_text

    def test_processes_defined(self, nf_text):
        expected_processes = [
            "PREPARE_EXOMEDEPTH_TARGETS",
            "BATCH_BAMS_EXOMEDEPTH",
            "GEN_EXOMEDEPTH_COUNTS",
            "MERGE_EXOMEDEPTH_COUNTS",
            "RUN_EXOMEDEPTH",
            "CONVERT_EXOMEDEPTH_TO_VCF",
        ]
        for proc in expected_processes:
            assert f"process {proc}" in nf_text, f"Module must define process {proc}"

    def test_process_labels(self, nf_text):
        assert "label 'bedtools'" in nf_text
        assert "label 'canoes_R'" in nf_text
        assert "label 'exomedepth'" in nf_text
        assert "label 'pysam'" in nf_text


# ===========================================================================
# 2. VCF Converter CLI & Formatting
# ===========================================================================

class TestExomeDepthConverter:

    def test_file_exists(self):
        assert os.path.isfile(CONVERTER), "bin/exomedepth_csv_to_vcf.py must exist"

    def test_has_main_and_cli_args(self, converter_text):
        assert "def main(" in converter_text
        assert "--input_file" in converter_text
        assert "--output_dir" in converter_text
        assert "--sample_id" in converter_text
        assert "--fai_file" in converter_text
        assert "--min_bf" in converter_text

    def test_vcf_header_contains_exomedepth_meta(self, converter_text):
        assert "##source=ExomeDepth" in converter_text
        assert "TOOL=EXOMEDEPTH" in converter_text
        assert "ID=BF" in converter_text
        assert "ID=READS_RATIO" in converter_text

    def test_conversion_execution(self):
        # Dynamically import exomedepth_csv_to_vcf
        spec = importlib.util.spec_from_file_location("exomedepth_csv_to_vcf", CONVERTER)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = os.path.join(tmpdir, "test.csv")
            with open(csv_path, "w") as f:
                f.write(
                    "chromosome,start,end,type,BF,reads.expected,reads.observed,reads.ratio,nexons\n"
                    "1,1000,2000,duplication,25.5,50.0,100.0,2.0,5\n"
                    "1,3000,4000,deletion,5.0,80.0,20.0,0.25,3\n"
                )

            fai_path = os.path.join(tmpdir, "ref.fa.fai")
            with open(fai_path, "w") as f:
                f.write("chr1\t248956422\t112\t70\t71\n")

            out_vcf = mod.convert_csv_to_vcf(
                input_csv=csv_path,
                sample_id="SAMPLE1",
                output_dir=tmpdir,
                fai_file=fai_path,
                min_bf=10.0,
            )

            assert os.path.isfile(out_vcf)
            with open(out_vcf) as fh:
                lines = fh.readlines()

            # Verify header lines
            assert any(l.startswith("##fileformat=VCFv4.2") for l in lines)
            assert any(l.startswith("##source=ExomeDepth") for l in lines)
            assert any(l.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1") for l in lines)

            # Variant records (non-header)
            var_lines = [l for l in lines if not l.startswith("#")]
            assert len(var_lines) == 2

            # Line 1: duplication with BF=25.5 >= 10.0 -> PASS
            cols1 = var_lines[0].split("\t")
            assert cols1[0] in ("1", "chr1")
            assert cols1[1] == "1000"
            assert cols1[4] == "<DUP>"
            assert cols1[6] == "PASS"
            assert "TOOL=EXOMEDEPTH" in cols1[7]
            assert "BF=25.5" in cols1[7]

            # Line 2: deletion with BF=5.0 < 10.0 -> LowQual
            cols2 = var_lines[1].split("\t")
            assert cols2[0] in ("1", "chr1")
            assert cols2[1] == "3000"
            assert cols2[4] == "<DEL>"
            assert cols2[6] == "LowQual"
            assert "BF=5.0" in cols2[7]


# ===========================================================================
# 3. Machine Learning & Feature Extraction Integration
# ===========================================================================

class TestExomeDepthMLIntegration:

    def test_feature_extraction_supported_callers(self):
        fe_path = os.path.join(REPO_ROOT, 'bin', 'feature_extraction.py')
        spec = importlib.util.spec_from_file_location("feature_extraction", fe_path)
        fe = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fe)
        assert 'exomedepth' in fe.SUPPORTED_CALLERS

    def test_train_xgboost_supported_callers(self):
        tx_path = os.path.join(REPO_ROOT, 'benchmarks', 'model_comparisons', 'train_xgboost.py')
        spec = importlib.util.spec_from_file_location("train_xgboost", tx_path)
        tx = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(tx)
        assert 'exomedepth' in tx.SUPPORTED_CALLERS
