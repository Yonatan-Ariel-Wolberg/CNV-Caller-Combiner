#!/usr/bin/env python3
"""
Tests for modules/ml/modules-feature-extraction.nf.

Validates the Nextflow module by reading the .nf file as text and asserting on
process/workflow content via string checks (consistent with the repo's existing
test practice).
"""

import os
import pytest

NF_PATH = os.path.join(
    os.path.dirname(__file__),
    '..',
    'modules',
    'ml',
    'modules-feature-extraction.nf',
)
MAIN_NF = os.path.join(os.path.dirname(__file__), '..', 'main.nf')


@pytest.fixture(scope='module')
def nf_text():
    with open(NF_PATH) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def main_nf_text():
    with open(MAIN_NF) as fh:
        return fh.read()


# ===========================================================================
# Module-level structure
# ===========================================================================

class TestModuleStructure:

    def test_file_exists(self):
        assert os.path.isfile(NF_PATH), (
            "modules/ml/modules-feature-extraction.nf must exist"
        )

    def test_dsl2_enabled(self, nf_text):
        assert 'nextflow.enable.dsl=2' in nf_text

    def test_extract_features_process(self, nf_text):
        assert 'process EXTRACT_FEATURES' in nf_text

    def test_feature_extraction_workflow(self, nf_text):
        assert 'workflow FEATURE_EXTRACTION' in nf_text

    def test_tag_uses_sample_id(self, nf_text):
        assert 'tag "${sample_id}"' in nf_text

    def test_label_feature_extraction(self, nf_text):
        assert "label 'feature_extraction'" in nf_text

    def test_publishdir_present(self, nf_text):
        assert 'publishDir' in nf_text


# ===========================================================================
# Process I/O
# ===========================================================================

class TestProcessIO:

    def test_input_sample_id(self, nf_text):
        assert 'sample_id' in nf_text

    def test_input_merged_vcf(self, nf_text):
        assert 'merged_vcf' in nf_text

    def test_input_tool_vcfs_str(self, nf_text):
        assert 'tool_vcfs_str' in nf_text

    def test_input_merger_mode(self, nf_text):
        assert 'merger_mode' in nf_text

    def test_input_bam_file(self, nf_text):
        assert 'bam_file' in nf_text

    def test_input_bed_file(self, nf_text):
        assert 'bed_file' in nf_text

    def test_input_mappability_file(self, nf_text):
        assert 'mappability_file' in nf_text

    def test_input_indelible_counts(self, nf_text):
        assert 'indelible_counts' in nf_text

    def test_output_features_tsv(self, nf_text):
        assert 'features_tsv' in nf_text and '_features.tsv' in nf_text

    def test_emit_features_tsv(self, nf_text):
        assert 'emit: features_tsv' in nf_text


# ===========================================================================
# Script block content
# ===========================================================================

class TestScriptBlock:

    def test_calls_feature_extraction_py(self, nf_text):
        assert 'feature_extraction.py' in nf_text

    def test_merged_vcf_arg(self, nf_text):
        assert '--merged_vcf' in nf_text

    def test_output_arg(self, nf_text):
        assert '--output' in nf_text

    def test_merger_mode_arg(self, nf_text):
        assert '--merger_mode' in nf_text

    def test_sample_id_arg(self, nf_text):
        assert '--sample_id' in nf_text

    def test_optional_args_conditional(self, nf_text):
        """Optional file arguments should be conditionally constructed."""
        assert 'bam_arg' in nf_text or '--bam_file' in nf_text
        assert 'bed_arg' in nf_text or '--bed_file' in nf_text


# ===========================================================================
# main.nf integration
# ===========================================================================

class TestMainNfIntegration:

    def test_includes_feature_extraction_module(self, main_nf_text):
        assert "include { FEATURE_EXTRACTION }" in main_nf_text or \
               "include { FEATURE_EXTRACTION}" in main_nf_text, (
            "main.nf must include FEATURE_EXTRACTION from modules-feature-extraction.nf"
        )

    def test_run_feature_extraction_workflow(self, main_nf_text):
        assert 'RUN_FEATURE_EXTRACTION' in main_nf_text

    def test_feature_extraction_case(self, main_nf_text):
        assert "case['feature_extraction']" in main_nf_text

    def test_feature_extraction_in_help_message(self, main_nf_text):
        assert 'feature_extraction' in main_nf_text

    def test_merged_vcf_dir_param_used(self, main_nf_text):
        assert 'merged_vcf_dir' in main_nf_text


# ===========================================================================
# nextflow.config integration
# ===========================================================================

class TestNextflowConfigIntegration:

    def test_feature_extraction_label_in_config(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'nextflow.config')
        with open(config_path) as fh:
            config_text = fh.read()
        assert "withLabel: 'feature_extraction'" in config_text, (
            "nextflow.config must define a 'feature_extraction' process label"
        )
