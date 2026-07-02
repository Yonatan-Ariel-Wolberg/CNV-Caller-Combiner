#!/usr/bin/env python3
"""
Tests for modules/ml/modules-train.nf and its integration into main.nf.

Validates:
  1. modules-train.nf exists and has the correct DSL2 / process / workflow structure.
  2. TRAIN_XGBOOST process inputs / outputs are correct.
  3. The script block calls train_xgboost.py with the correct arguments.
  4. main.nf includes the TRAIN module, defines RUN_TRAIN, and has a
     case['train'] block that reads --features_dir and --truth_labels params.
  5. nextflow.config defines a 'train' process label.
  6. params/general/params-train.json exists and contains valid JSON with the correct keys.
  7. bin/train_xgboost.py exposes a main() function and uses argparse.
"""

import json
import os
import re

import pytest

REPO_ROOT   = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
NF_PATH     = os.path.join(REPO_ROOT, 'modules', 'ml', 'modules-train.nf')
MAIN_NF     = os.path.join(REPO_ROOT, 'main.nf')
CONFIG_PATH = os.path.join(REPO_ROOT, 'nextflow.config')
PARAMS_FILE = os.path.join(REPO_ROOT, 'params', 'general', 'params-train.json')
SCRIPT_PATH = os.path.join(REPO_ROOT, 'bin', 'train_xgboost.py')
TRAIN_DEF = os.path.join(REPO_ROOT, 'bin', 'train.def')


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def nf_text():
    with open(NF_PATH) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def main_text():
    with open(MAIN_NF) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def config_text():
    with open(CONFIG_PATH) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def script_text():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


# ===========================================================================
# 1. Module structure
# ===========================================================================

class TestModuleStructure:
    """modules-train.nf must be a well-formed DSL2 module."""

    def test_file_exists(self):
        assert os.path.isfile(NF_PATH), (
            'modules/ml/modules-train.nf must exist'
        )

    def test_dsl2_enabled(self, nf_text):
        assert 'nextflow.enable.dsl=2' in nf_text

    def test_has_train_xgboost_process(self, nf_text):
        assert 'process TRAIN_XGBOOST' in nf_text, (
            'Module must define process TRAIN_XGBOOST'
        )

    def test_has_train_workflow(self, nf_text):
        assert 'workflow TRAIN' in nf_text, (
            'Module must define workflow TRAIN'
        )

    def test_uses_train_label(self, nf_text):
        assert "label 'train'" in nf_text, (
            "Process must use label 'train'"
        )

    def test_publishdir_present(self, nf_text):
        assert 'publishDir' in nf_text

    def test_publishdir_out_train(self, nf_text):
        assert 'out_TRAIN' in nf_text


# ===========================================================================
# 2. Process inputs / outputs
# ===========================================================================

class TestProcessIO:
    """TRAIN_XGBOOST must accept and emit the correct channels."""

    def test_input_features_tsv_files(self, nf_text):
        assert 'features_tsv_files' in nf_text

    def test_input_truth_labels(self, nf_text):
        assert 'truth_labels' in nf_text

    def test_input_optional_probes_bed(self, nf_text):
        assert 'probes_bed' in nf_text

    def test_output_model(self, nf_text):
        assert 'cnv_model.json' in nf_text

    def test_output_report(self, nf_text):
        assert 'training_report.txt' in nf_text

    def test_output_roc_plot(self, nf_text):
        assert 'roc_curve.svg' in nf_text

    def test_output_pr_plot(self, nf_text):
        assert 'pr_curve.svg' in nf_text

    def test_output_shap_values(self, nf_text):
        assert 'shap_values.tsv' in nf_text

    def test_emit_model(self, nf_text):
        assert 'emit: model' in nf_text

    def test_emit_report(self, nf_text):
        assert 'emit: report' in nf_text


# ===========================================================================
# 3. Script block content
# ===========================================================================

class TestScriptBlock:
    """The script block must call train_xgboost.py with the expected arguments."""

    def test_calls_train_script(self, nf_text):
        assert 'train_xgboost.py' in nf_text

    def test_passes_features_dir(self, nf_text):
        assert '--features_dir' in nf_text

    def test_passes_truth_labels(self, nf_text):
        assert '--truth_labels' in nf_text

    def test_supports_optional_probes_bed(self, nf_text):
        assert '--probes_bed' in nf_text

    def test_passes_output_model(self, nf_text):
        assert '--output_model' in nf_text

    def test_passes_output_report(self, nf_text):
        assert '--output_report' in nf_text

    def test_passes_roc_and_pr_outputs(self, nf_text):
        assert '--output_roc_plot' in nf_text
        assert '--output_pr_plot' in nf_text
        assert '--output_roc_data' in nf_text
        assert '--output_pr_data' in nf_text

    def test_passes_shap_outputs(self, nf_text):
        assert '--output_shap_values' in nf_text
        assert '--output_shap_summary_plot' in nf_text
        assert '--output_shap_beeswarm_plot' in nf_text


# ===========================================================================
# 4. Workflow TRAIN take/main/emit structure
# ===========================================================================

class TestWorkflowTrain:
    """The TRAIN workflow must expose the expected take/main/emit structure."""

    def _train_block(self, nf_text):
        m = re.search(r'workflow TRAIN\s*\{(.+)', nf_text, re.DOTALL)
        assert m, 'TRAIN workflow block not found'
        return m.group(1)

    def test_workflow_has_take(self, nf_text):
        block = self._train_block(nf_text)
        assert 'take:' in block

    def test_workflow_has_main(self, nf_text):
        block = self._train_block(nf_text)
        assert 'main:' in block

    def test_workflow_calls_process(self, nf_text):
        block = self._train_block(nf_text)
        assert 'TRAIN_XGBOOST' in block
        assert 'probes_bed_ch' in block

    def test_workflow_has_emit(self, nf_text):
        block = self._train_block(nf_text)
        assert 'emit:' in block

    def test_workflow_emits_model(self, nf_text):
        block = self._train_block(nf_text)
        assert 'model' in block

    def test_workflow_emits_report(self, nf_text):
        block = self._train_block(nf_text)
        assert 'report' in block

    def test_workflow_emits_new_artifacts(self, nf_text):
        block = self._train_block(nf_text)
        assert 'roc_plot' in block
        assert 'pr_plot' in block
        assert 'shap_values' in block
        assert 'shap_summary_bar_plot' in block
        assert 'shap_beeswarm_plot' in block

    def test_workflow_collects_features(self, nf_text):
        """Feature TSVs must be collected so the process sees all samples."""
        block = self._train_block(nf_text)
        assert '.collect()' in block


# ===========================================================================
# 5. main.nf integration
# ===========================================================================

class TestMainNfIntegration:
    """main.nf must include and wire up the TRAIN module correctly."""

    def test_includes_train_module(self, main_text):
        assert "include { TRAIN }" in main_text or \
               "include { TRAIN}" in main_text, (
            "main.nf must include { TRAIN } from modules-train.nf"
        )
        assert 'modules-train.nf' in main_text, (
            "main.nf include statement must reference modules-train.nf"
        )

    def test_run_train_workflow_defined(self, main_text):
        assert 'workflow RUN_TRAIN' in main_text, (
            "main.nf must define 'workflow RUN_TRAIN'"
        )

    def test_train_case_in_switch(self, main_text):
        assert "case['train']" in main_text, (
            "main.nf switch block must contain case['train']"
        )

    def _train_case_block(self, main_text):
        m = re.search(r"case\['train'\](.+?)break", main_text, re.DOTALL)
        assert m, "case['train'] block not found"
        return m.group(1)

    def test_train_case_reads_features_dir(self, main_text):
        block = self._train_case_block(main_text)
        assert 'features_dir' in block, (
            "case['train'] must read --features_dir param"
        )

    def test_train_case_reads_truth_labels(self, main_text):
        block = self._train_case_block(main_text)
        assert 'truth_labels' in block, (
            "case['train'] must read --truth_labels param"
        )

    def test_train_case_reads_optional_probes_bed(self, main_text):
        block = self._train_case_block(main_text)
        assert 'probes_bed' in block

    def test_train_case_calls_run_train(self, main_text):
        block = self._train_case_block(main_text)
        assert 'RUN_TRAIN' in block, (
            "case['train'] must call RUN_TRAIN"
        )
        assert 'ch_probes' in block

    def test_train_in_help_message(self, main_text):
        assert '--workflow train' in main_text, (
            "The default error message must list '--workflow train' as an option"
        )


# ===========================================================================
# 6. nextflow.config integration
# ===========================================================================

class TestNextflowConfigIntegration:
    """nextflow.config must define a 'train' process label with a container."""

    def test_train_label_in_config(self, config_text):
        assert "withLabel: 'train'" in config_text, (
            "nextflow.config must define a 'train' process label"
        )

    def test_train_label_has_container(self, config_text):
        # Extract the train label block and verify it defines a container
        m = re.search(
            r"withLabel:\s*'train'\s*\{([^}]+)\}",
            config_text,
            re.DOTALL,
        )
        assert m, "withLabel: 'train' block not found in nextflow.config"
        block = m.group(1)
        assert 'container' in block, (
            "withLabel: 'train' block must define a container"
        )

    def test_train_label_has_memory(self, config_text):
        m = re.search(
            r"withLabel:\s*'train'\s*\{([^}]+)\}",
            config_text,
            re.DOTALL,
        )
        assert m, "withLabel: 'train' block not found"
        assert 'memory' in m.group(1)

    def test_train_label_has_cpus(self, config_text):
        m = re.search(
            r"withLabel:\s*'train'\s*\{([^}]+)\}",
            config_text,
            re.DOTALL,
        )
        assert m, "withLabel: 'train' block not found"
        assert 'cpus' in m.group(1)

    def test_train_label_uses_xgboost_biocontainer(self, config_text):
        m = re.search(
            r"withLabel:\s*'train'\s*\{([^}]+)\}",
            config_text,
            re.DOTALL,
        )
        assert m, "withLabel: 'train' block not found"
        assert "docker://quay.io/condaforge/mambaforge:24.9.2-0" in m.group(1)

    def test_train_label_does_not_use_legacy_py27_xgboost_container(self, config_text):
        m = re.search(
            r"withLabel:\s*'train'\s*\{([^}]+)\}",
            config_text,
            re.DOTALL,
        )
        assert m, "withLabel: 'train' block not found"
        container_line = re.search(r"container\s*=\s*'([^']+)'", m.group(1))
        assert container_line, "train label must define container assignment"
        assert "xgboost:0.6a2--py27_0" not in container_line.group(1)


# ===========================================================================
# 7. params/general/params-train.json
# ===========================================================================

class TestParamsFile:
    """params/general/params-train.json must exist and be valid JSON with correct keys."""

    def test_params_file_exists(self):
        assert os.path.isfile(PARAMS_FILE), (
            'params/general/params-train.json must exist'
        )

    def test_params_file_valid_json(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert isinstance(data, dict), 'params-train.json must be a JSON object'

    def test_params_file_workflow_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert data.get('workflow') == 'train', (
            'params-train.json must contain "workflow": "train"'
        )

    def test_params_file_features_dir_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert 'features_dir' in data, (
            'params-train.json must include a features_dir key'
        )

    def test_params_file_truth_labels_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert 'truth_labels' in data, (
            'params-train.json must include a truth_labels key'
        )


# ===========================================================================
# 8. Apptainer recipe for train dependencies
# ===========================================================================

class TestTrainApptainerRecipe:
    """bin/train.def must define the pinned train dependencies."""

    def test_train_def_exists(self):
        assert os.path.isfile(TRAIN_DEF), 'bin/train.def must exist'

    def test_train_def_uses_expected_base_image(self):
        with open(TRAIN_DEF) as fh:
            text = fh.read()
        assert 'Bootstrap: docker' in text
        assert 'From: quay.io/condaforge/mambaforge:24.9.2-0' in text

    def test_train_def_pins_expected_conda_packages(self):
        with open(TRAIN_DEF) as fh:
            text = fh.read()
        assert 'xgboost=2.1.4' in text
        assert 'scikit-learn=1.6.1' in text
        assert 'imbalanced-learn=0.13.0' in text
        assert 'pandas=2.2.3' in text
        assert 'numpy=2.2.3' in text


# ===========================================================================
# 8. bin/train_xgboost.py – CLI integration
# ===========================================================================

class TestTrainScriptCLI:
    """bin/train_xgboost.py must be a proper CLI script."""

    def test_has_argparse_import(self, script_text):
        assert 'import argparse' in script_text

    def test_has_main_function(self, script_text):
        assert 'def main(' in script_text

    def test_has_name_main_guard(self, script_text):
        assert "if __name__ == '__main__'" in script_text or \
               'if __name__ == "__main__"' in script_text, (
            "bin/train_xgboost.py must have an if __name__ == '__main__' guard"
        )

    def test_main_accepts_features_dir(self, script_text):
        assert '--features_dir' in script_text

    def test_main_accepts_truth_labels(self, script_text):
        assert '--truth_labels' in script_text

    def test_main_accepts_output_model(self, script_text):
        assert '--output_model' in script_text

    def test_main_accepts_output_report(self, script_text):
        assert '--output_report' in script_text

    def test_main_saves_model(self, script_text):
        assert 'save_model' in script_text

    def test_main_writes_report(self, script_text):
        assert 'output_report' in script_text

    def test_main_loads_feature_tsvs(self, script_text):
        """Script must discover *_features.tsv files from features_dir."""
        assert '_features.tsv' in script_text

    def test_main_merges_truth_labels(self, script_text):
        """Script must join features with truth labels."""
        assert 'truth_label' in script_text


# ===========================================================================
# 9. README training documentation
# ===========================================================================

class TestReadmeTrainTruthLabelsDocumentation:
    """README should document required truth-label fields for training."""

    README_PATH = os.path.join(REPO_ROOT, 'README.md')

    @pytest.fixture(scope='class')
    def readme_text(self):
        with open(self.README_PATH) as fh:
            return fh.read()

    def test_readme_mentions_truth_label_tsv_requirements(self, readme_text):
        assert 'Truth-label TSV requirements' in readme_text

    def test_readme_lists_required_truth_label_columns(self, readme_text):
        assert '`sample_id`' in readme_text
        assert '`chrom`' in readme_text
        assert '`start`' in readme_text
        assert '`end`' in readme_text
        assert '`cnv_type`' in readme_text
        assert '`truth_label`' in readme_text

    def test_readme_documents_truth_label_value(self, readme_text):
        section_start = readme_text.find('Truth-label TSV requirements')
        assert section_start != -1
        heading_match = re.search(r'\n#{1,6}(\s|\n|$)', readme_text[section_start + 1:])
        section_end = (
            section_start + 1 + heading_match.start()
            if heading_match else len(readme_text)
        )
        section = readme_text[section_start:section_end]
        assert 'Example:' in section
        assert re.search(r'truth_label\s*=\s*1', section)
        assert re.search(r'truth_label\s*=\s*0', section)
