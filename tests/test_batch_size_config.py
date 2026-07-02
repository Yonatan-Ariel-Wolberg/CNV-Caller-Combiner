#!/usr/bin/env python3
"""
Tests for batch-size configuration of CANOES and XHMM modules.

Validates that:
  1. The CANOES module default batch size is 100 (not 200), matching the
     medium-profile value and README documentation.
  2. The XHMM module default batch size is 50, matching the medium-profile
     value and README documentation.
  3. nextflow.config declares global defaults of 100 / 50 for these params,
     ensuring a sensible value is always available even when no cohort-size
     profile is selected.
  4. The cohort-size profiles in nextflow.config set the documented
     per-profile values for both parameters.
  5. The example params-JSON files carry the correct batch-size values.
"""

import os
import re
import json

# Repository root relative to this test file
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

CANOES_MODULE   = os.path.join(REPO_ROOT, 'modules', 'callers', 'modules-canoes.nf')
XHMM_MODULE     = os.path.join(REPO_ROOT, 'modules', 'callers', 'modules-xhmm.nf')
NEXTFLOW_CONFIG = os.path.join(REPO_ROOT, 'nextflow.config')
PARAMS_CANOES   = os.path.join(REPO_ROOT, 'params', 'general', 'params-canoes.json')
PARAMS_XHMM     = os.path.join(REPO_ROOT, 'params', 'general', 'params-xhmm.json')
XHMM_CONF       = os.path.join(REPO_ROOT, 'params', 'xhmm.conf')


# ---------------------------------------------------------------------------
# Module-file defaults
# ---------------------------------------------------------------------------

class TestModuleDefaults:
    """The module-level params.X = value fallbacks must match the documented
    medium-profile defaults."""

    def test_canoes_module_default_is_100(self):
        """params.canoes_batch_size in modules-canoes.nf must default to 100."""
        with open(CANOES_MODULE) as f:
            content = f.read()
        # Match: params.canoes_batch_size = <NUMBER>
        match = re.search(r'params\.canoes_batch_size\s*=\s*(\d+)', content)
        assert match is not None, \
            "params.canoes_batch_size not found in modules-canoes.nf"
        value = int(match.group(1))
        assert value == 100, (
            f"Expected params.canoes_batch_size = 100 in modules-canoes.nf "
            f"(medium-profile default), got {value}"
        )

    def test_xhmm_module_default_is_50(self):
        """params.xhmm_batch_size in modules-xhmm.nf must default to 50."""
        with open(XHMM_MODULE) as f:
            content = f.read()
        match = re.search(r'params\.xhmm_batch_size\s*=\s*(\d+)', content)
        assert match is not None, \
            "params.xhmm_batch_size not found in modules-xhmm.nf"
        value = int(match.group(1))
        assert value == 50, (
            f"Expected params.xhmm_batch_size = 50 in modules-xhmm.nf "
            f"(medium-profile default), got {value}"
        )


# ---------------------------------------------------------------------------
# nextflow.config global-params defaults
# ---------------------------------------------------------------------------

class TestNextflowConfigGlobalDefaults:
    """The main params { } block must declare sensible global defaults so users
    who select no cohort-size profile still get a reasonable value."""

    def _load_config(self):
        with open(NEXTFLOW_CONFIG) as f:
            return f.read()

    def test_global_canoes_batch_size_is_100(self):
        """Global params block should declare canoes_batch_size = 100."""
        content = self._load_config()
        # The global params block appears before the 'profiles {' block.
        global_section = content.split('profiles {')[0]
        match = re.search(r'canoes_batch_size\s*=\s*(\d+)', global_section)
        assert match is not None, \
            "canoes_batch_size not found in the global params block of nextflow.config"
        value = int(match.group(1))
        assert value == 100, (
            f"Global canoes_batch_size should be 100 (medium default), got {value}"
        )

    def test_global_xhmm_batch_size_is_50(self):
        """Global params block should declare xhmm_batch_size = 50."""
        content = self._load_config()
        global_section = content.split('profiles {')[0]
        match = re.search(r'xhmm_batch_size\s*=\s*(\d+)', global_section)
        assert match is not None, \
            "xhmm_batch_size not found in the global params block of nextflow.config"
        value = int(match.group(1))
        assert value == 50, (
            f"Global xhmm_batch_size should be 50 (medium default), got {value}"
        )

    def test_global_help_default_is_false(self):
        """Global params block should declare help = false."""
        content = self._load_config()
        global_section = content.split('profiles {')[0]
        match = re.search(r'help\s*=\s*(\S+)', global_section)
        assert match is not None, \
            "help not found in the global params block of nextflow.config"
        assert match.group(1) == 'false', (
            f"Global help should be false by default, got {match.group(1)}"
        )

    def test_global_test_size_default_is_minus_one(self):
        """Global params block should declare test_size = -1."""
        content = self._load_config()
        global_section = content.split('profiles {')[0]
        match = re.search(r'test_size\s*=\s*(-?\d+)', global_section)
        assert match is not None, \
            "test_size not found in the global params block of nextflow.config"
        assert int(match.group(1)) == -1, (
            f"Global test_size should be -1 by default, got {match.group(1)}"
        )


# ---------------------------------------------------------------------------
# nextflow.config cohort-size profile values
# ---------------------------------------------------------------------------

class TestNextflowConfigProfiles:
    """The three cohort-size profiles must set the documented batch-size values
    so that the correct value is used when a profile is selected."""

    def _load_config(self):
        with open(NEXTFLOW_CONFIG) as f:
            return f.read()

    @staticmethod
    def _profile_block(content, profile_name):
        """Extract the text of a named profile block (including nested braces)."""
        pattern = rf'{re.escape(profile_name)}\s*\{{'
        m = re.search(pattern, content)
        assert m is not None, f"Profile '{profile_name}' not found in nextflow.config"
        start = m.end() - 1          # position of opening '{'
        depth = 0
        for i, ch in enumerate(content[start:], start):
            if ch == '{':
                depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0:
                    return content[start: i + 1]
        raise AssertionError(f"Unmatched braces for profile '{profile_name}'")

    @staticmethod
    def _param_value(block, param_name):
        m = re.search(rf'{re.escape(param_name)}\s*=\s*(\d+)', block)
        assert m is not None, \
            f"'{param_name}' not found in profile block:\n{block}"
        return int(m.group(1))

    # --- small profile ---

    def test_small_profile_canoes_batch_size(self):
        block = self._profile_block(self._load_config(), 'small')
        assert self._param_value(block, 'canoes_batch_size') == 50

    def test_small_profile_xhmm_batch_size(self):
        block = self._profile_block(self._load_config(), 'small')
        assert self._param_value(block, 'xhmm_batch_size') == 25

    # --- medium profile ---

    def test_medium_profile_canoes_batch_size(self):
        block = self._profile_block(self._load_config(), 'medium')
        assert self._param_value(block, 'canoes_batch_size') == 100

    def test_medium_profile_xhmm_batch_size(self):
        block = self._profile_block(self._load_config(), 'medium')
        assert self._param_value(block, 'xhmm_batch_size') == 50

    # --- large profile ---

    def test_large_profile_canoes_batch_size(self):
        block = self._profile_block(self._load_config(), 'large')
        assert self._param_value(block, 'canoes_batch_size') == 200

    def test_large_profile_xhmm_batch_size(self):
        block = self._profile_block(self._load_config(), 'large')
        assert self._param_value(block, 'xhmm_batch_size') == 100

    # --- ordering: small < medium < large ---

    def test_profile_canoes_batch_sizes_are_ordered(self):
        content = self._load_config()
        small  = self._param_value(self._profile_block(content, 'small'),  'canoes_batch_size')
        medium = self._param_value(self._profile_block(content, 'medium'), 'canoes_batch_size')
        large  = self._param_value(self._profile_block(content, 'large'),  'canoes_batch_size')
        assert small < medium < large, (
            f"Expected small < medium < large for canoes_batch_size, "
            f"got {small} / {medium} / {large}"
        )

    def test_profile_xhmm_batch_sizes_are_ordered(self):
        content = self._load_config()
        small  = self._param_value(self._profile_block(content, 'small'),  'xhmm_batch_size')
        medium = self._param_value(self._profile_block(content, 'medium'), 'xhmm_batch_size')
        large  = self._param_value(self._profile_block(content, 'large'),  'xhmm_batch_size')
        assert small < medium < large, (
            f"Expected small < medium < large for xhmm_batch_size, "
            f"got {small} / {medium} / {large}"
        )


# ---------------------------------------------------------------------------
# Example params-JSON files
# ---------------------------------------------------------------------------

class TestParamsJsonFiles:
    """The example params-file templates must specify reasonable batch sizes."""

    def test_params_canoes_json_batch_size(self):
        with open(PARAMS_CANOES) as f:
            data = json.load(f)
        assert 'canoes_batch_size' in data, \
            "params-canoes.json is missing canoes_batch_size"
        assert data['canoes_batch_size'] == 100, (
            f"params-canoes.json: expected canoes_batch_size=100, "
            f"got {data['canoes_batch_size']}"
        )

    def test_params_xhmm_json_batch_size(self):
        with open(PARAMS_XHMM) as f:
            data = json.load(f)
        assert 'xhmm_batch_size' in data, \
            "params-xhmm.json is missing xhmm_batch_size"
        assert data['xhmm_batch_size'] == 50, (
            f"params-xhmm.json: expected xhmm_batch_size=50, "
            f"got {data['xhmm_batch_size']}"
        )

    def test_params_xhmm_json_uses_repo_xhmm_conf(self):
        with open(PARAMS_XHMM) as f:
            data = json.load(f)
        assert data.get('xhmm_conf') == "params/xhmm.conf", (
            "params-xhmm.json should reference the repository template "
            "config at params/xhmm.conf"
        )

    def test_repo_xhmm_conf_exists_and_has_nine_fields(self):
        assert os.path.isfile(XHMM_CONF), "params/xhmm.conf must exist"
        with open(XHMM_CONF) as f:
            line = f.read().strip()
        fields = line.split()
        assert len(fields) == 9, (
            "params/xhmm.conf must contain one whitespace-delimited line "
            "with 9 XHMM HMM parameters"
        )


# ---------------------------------------------------------------------------
# Normalisation-design guarantees (structure checks)
# ---------------------------------------------------------------------------

class TestNormalisationDesign:
    """Verify that the pipeline merges all batches BEFORE normalisation for
    each caller that uses batching (CANOES and XHMM).

    These are structural smoke tests that check that the merge/combine step
    precedes the normalisation step in the workflow source code."""

    def test_canoes_merge_precedes_run_canoes(self):
        """MERGE_READ_COUNTS must appear before RUN_CANOES in the CANOES
        workflow block so that all batches are combined before normalisation."""
        with open(CANOES_MODULE) as f:
            content = f.read()
        # Find the CANOES workflow block
        wf_match = re.search(r'workflow CANOES \{', content)
        assert wf_match, "workflow CANOES not found in modules-canoes.nf"
        wf_block = content[wf_match.start():]
        merge_pos = wf_block.find('MERGE_READ_COUNTS')
        run_pos   = wf_block.find('RUN_CANOES')
        assert merge_pos != -1, \
            "MERGE_READ_COUNTS not found in workflow CANOES"
        assert run_pos != -1, \
            "RUN_CANOES not found in workflow CANOES"
        assert merge_pos < run_pos, (
            "MERGE_READ_COUNTS must appear before RUN_CANOES in the CANOES "
            "workflow to ensure all batches are merged before normalisation"
        )

    def test_canoes_sample_list_channel_has_ifempty_guard(self):
        """sample_list_ch derivation must fail fast if merged outputs are empty."""
        with open(CANOES_MODULE) as f:
            content = f.read()
        assert "sample_list_ch = MERGE_READ_COUNTS.out.chr_reads_cov" in content, (
            "Expected sample_list_ch assignment in CANOES workflow"
        )
        sample_assign_match = re.search(
            r"sample_list_ch\s*=\s*MERGE_READ_COUNTS\.out\.chr_reads_cov(.+?)\.map\s*\{",
            content,
            re.DOTALL,
        )
        assert sample_assign_match is not None, (
            "Could not locate sample_list_ch derivation chain in CANOES workflow"
        )
        chain_text = sample_assign_match.group(1)
        assert ".ifEmpty" in chain_text, (
            "sample_list_ch derivation should include .ifEmpty to avoid hanging "
            "or unclear failures when MERGE_READ_COUNTS emits nothing."
        )

    def test_xhmm_combine_precedes_filter_samples(self):
        """COMBINE_DOC must appear before FILTER_SAMPLES (which starts the
        normalisation pipeline) in the XHMM workflow block."""
        with open(XHMM_MODULE) as f:
            content = f.read()
        wf_match = re.search(r'workflow XHMM \{', content)
        assert wf_match, "workflow XHMM not found in modules-xhmm.nf"
        wf_block = content[wf_match.start():]
        combine_pos = wf_block.find('COMBINE_DOC')
        filter_pos  = wf_block.find('FILTER_SAMPLES')
        assert combine_pos != -1, \
            "COMBINE_DOC not found in workflow XHMM"
        assert filter_pos != -1, \
            "FILTER_SAMPLES not found in workflow XHMM"
        assert combine_pos < filter_pos, (
            "COMBINE_DOC must appear before FILTER_SAMPLES in the XHMM "
            "workflow to ensure all batches are merged before normalisation"
        )
