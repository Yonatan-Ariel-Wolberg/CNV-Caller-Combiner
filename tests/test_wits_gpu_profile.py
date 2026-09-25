#!/usr/bin/env python3
"""
Tests for GPU-node execution on the ZA-Wits-Core HPC cluster.

Validates that nextflow.config is correctly configured for:
  1. A `wits_gpu` profile exists in the profiles block.
  2. The profile configures the `RUN_PARABRICKS` process with the correct
     SLURM GPU resource request (`--gres=gpu:1 --reservation=gpu`) and the
     official NVIDIA Clara Parabricks container image.
  3. No GPU clusterOptions (--gres=gpu / --reservation=gpu) appear in the
     default process block or in the CPU-only `wits` profile—GPU resources
     are confined to the dedicated `wits_gpu` profile.
  4. The `wits_gpu` profile is compatible with the `wits` profile (i.e. can
     be combined on the command line as `-profile wits,wits_gpu`).
"""

import os
import re

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
NEXTFLOW_CONFIG = os.path.join(REPO_ROOT, 'nextflow.config')

PARABRICKS_CONTAINER = 'nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1'


def _read_config():
    with open(NEXTFLOW_CONFIG) as fh:
        return fh.read()


def _extract_block(content, start_pattern):
    """Return the body of the first block matching start_pattern { ... }.

    Uses brace counting to handle nested blocks correctly.
    """
    match = re.search(start_pattern + r'\s*\{', content)
    assert match, f"Block matching '{start_pattern}' not found in nextflow.config"
    start = match.end()  # position just after the opening {
    depth = 1
    pos = start
    while pos < len(content) and depth > 0:
        if content[pos] == '{':
            depth += 1
        elif content[pos] == '}':
            depth -= 1
        pos += 1
    return content[start:pos - 1]


def _profiles_section(content):
    """Return the body of the profiles { } block."""
    return _extract_block(content, r'\bprofiles')


def _wits_gpu_profile_section(content):
    """Return the text of the wits_gpu { ... } profile block."""
    profiles_text = _profiles_section(content)
    return _extract_block(profiles_text, r'\bwits_gpu')


def _wits_profile_section(content):
    """Return the text of the wits { ... } profile block (CPU-only)."""
    profiles_idx = content.find('profiles {')
    assert profiles_idx >= 0, "profiles { } block not found in nextflow.config"
    return _extract_block(content[profiles_idx:], r'\bwits\b')


def _global_process_section(content):
    """Return the text of the global process { } block (before profiles {)."""
    global_part = content.split('profiles {')[0]
    match = re.search(r'\bprocess\s*\{(.+)', global_part, re.DOTALL)
    assert match, "process block not found before profiles section"
    return match.group(1)


def _run_parabricks_block(profile_text):
    """Return the text of the withName: 'RUN_PARABRICKS' { } block."""
    return _extract_block(profile_text, r"withName:\s*'RUN_PARABRICKS'")


# ===========================================================================
# 1. wits_gpu profile exists
# ===========================================================================

class TestWitsGpuProfileExists:
    """The wits_gpu profile must be declared in nextflow.config."""

    def test_wits_gpu_profile_present(self):
        """profiles { } block must contain a wits_gpu { } section."""
        content = _read_config()
        profiles = _profiles_section(content)
        assert re.search(r'\bwits_gpu\s*\{', profiles), (
            "nextflow.config must define a 'wits_gpu' profile so that users can "
            "enable GPU-accelerated Parabricks execution on the NVIDIA L4 nodes "
            "of the ZA-Wits-Core cluster: -profile wits,wits_gpu"
        )

    def test_wits_gpu_profile_is_separate_from_wits(self):
        """wits_gpu must be a distinct profile, not merged into wits."""
        content = _read_config()
        profiles = _profiles_section(content)
        # Both profiles must appear as separate named blocks.
        wits_matches = list(re.finditer(r'\bwits\b\s*\{', profiles))
        wits_gpu_matches = list(re.finditer(r'\bwits_gpu\s*\{', profiles))
        assert len(wits_gpu_matches) >= 1, (
            "wits_gpu must be defined as its own profile block"
        )
        assert len(wits_matches) >= 1, (
            "wits (CPU-only) profile must still exist as a separate block"
        )


# ===========================================================================
# 2. RUN_PARABRICKS process configuration inside wits_gpu
# ===========================================================================

class TestRunParabricksProcess:
    """wits_gpu profile must configure RUN_PARABRICKS with GPU-specific settings."""

    def test_run_parabricks_block_present(self):
        """wits_gpu profile must contain a withName: 'RUN_PARABRICKS' block."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        assert re.search(r"withName:\s*'RUN_PARABRICKS'", gpu_profile), (
            "wits_gpu profile must configure the RUN_PARABRICKS process via "
            "withName: 'RUN_PARABRICKS' { ... }"
        )

    def test_run_parabricks_requests_gpu_gres(self):
        """RUN_PARABRICKS clusterOptions must request --gres=gpu:1."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        pb_block = _run_parabricks_block(gpu_profile)
        assert '--gres=gpu:1' in pb_block, (
            "RUN_PARABRICKS clusterOptions must include '--gres=gpu:1' so that "
            "SLURM allocates one GPU (NVIDIA L4 / Tesla K20Xm) for the job on "
            "the ZA-Wits-Core cluster"
        )

    def test_run_parabricks_requests_gpu_reservation(self):
        """RUN_PARABRICKS clusterOptions must include --reservation=gpu."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        pb_block = _run_parabricks_block(gpu_profile)
        assert '--reservation=gpu' in pb_block, (
            "RUN_PARABRICKS clusterOptions must include '--reservation=gpu' to "
            "target the dedicated GPU reservation on the Wits cluster"
        )

    def test_run_parabricks_container_is_parabricks(self):
        """RUN_PARABRICKS must use the official NVIDIA Clara Parabricks container."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        pb_block = _run_parabricks_block(gpu_profile)
        assert PARABRICKS_CONTAINER in pb_block, (
            f"RUN_PARABRICKS container must be '{PARABRICKS_CONTAINER}' — "
            "the official NVIDIA Clara Parabricks image that provides GPU-accelerated "
            "GATK-compatible alignment and variant calling"
        )

    def test_run_parabricks_cluster_options_line_present(self):
        """RUN_PARABRICKS block must contain a clusterOptions line."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        pb_block = _run_parabricks_block(gpu_profile)
        assert 'clusterOptions' in pb_block, (
            "RUN_PARABRICKS block must set clusterOptions to request GPU resources "
            "via the SLURM --gres mechanism"
        )

    def test_run_parabricks_both_gres_and_reservation_in_same_cluster_options(self):
        """Both --gres=gpu:1 and --reservation=gpu must be in a single clusterOptions."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        pb_block = _run_parabricks_block(gpu_profile)
        cluster_opts_line = next(
            (l for l in pb_block.splitlines() if 'clusterOptions' in l), None
        )
        assert cluster_opts_line is not None, "clusterOptions line not found in RUN_PARABRICKS"
        assert '--gres=gpu:1' in cluster_opts_line, (
            "--gres=gpu:1 must appear in the clusterOptions line"
        )
        assert '--reservation=gpu' in cluster_opts_line, (
            "--reservation=gpu must appear in the clusterOptions line"
        )


# ===========================================================================
# 3. GPU options must NOT appear in the default process block or wits profile
# ===========================================================================

class TestGpuOptionsIsolatedToGpuProfile:
    """GPU SLURM options must be confined to the wits_gpu profile."""

    def test_no_gres_gpu_in_default_process(self):
        """--gres=gpu must not appear in the global process block."""
        content = _read_config()
        proc = _global_process_section(content)
        assert '--gres=gpu' not in proc, (
            "The default process block must not contain '--gres=gpu'. "
            "GPU resources should only be requested in the wits_gpu profile."
        )

    def test_no_reservation_gpu_in_default_process(self):
        """--reservation=gpu must not appear in the global process block."""
        content = _read_config()
        proc = _global_process_section(content)
        assert '--reservation=gpu' not in proc, (
            "The default process block must not contain '--reservation=gpu'. "
            "GPU reservations should only be requested in the wits_gpu profile."
        )

    def test_no_gres_gpu_in_wits_cpu_profile(self):
        """--gres=gpu must not appear in the CPU-only wits profile."""
        content = _read_config()
        wits = _wits_profile_section(content)
        # Strip out any nested wits_gpu block before checking – it follows `wits`.
        # The _wits_profile_section helper already extracts only the `wits` block.
        assert '--gres=gpu' not in wits, (
            "The wits (CPU-only) profile must not contain '--gres=gpu'. "
            "GPU resources are confined to the wits_gpu profile."
        )

    def test_no_reservation_gpu_in_wits_cpu_profile(self):
        """--reservation=gpu must not appear in the CPU-only wits profile."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert '--reservation=gpu' not in wits, (
            "The wits (CPU-only) profile must not contain '--reservation=gpu'. "
            "GPU reservations are confined to the wits_gpu profile."
        )


# ===========================================================================
# 4. wits_gpu profile is composable with wits (does not duplicate executor setup)
# ===========================================================================

class TestWitsGpuComposability:
    """wits_gpu must not re-declare executor or SLURM settings already in wits."""

    def test_wits_gpu_does_not_set_executor(self):
        """wits_gpu must not set executor.name (already set by wits)."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        assert 'executor.name' not in gpu_profile, (
            "wits_gpu must not re-declare executor.name = 'slurm'. "
            "The executor is configured by the wits profile; wits_gpu is additive."
        )

    def test_wits_gpu_does_not_set_queue(self):
        """wits_gpu must not set process.queue (already set by wits)."""
        content = _read_config()
        gpu_profile = _wits_gpu_profile_section(content)
        assert 'process.queue' not in gpu_profile, (
            "wits_gpu must not re-declare process.queue. "
            "The queue is set by the wits profile; wits_gpu only overrides "
            "the RUN_PARABRICKS process."
        )
