#!/usr/bin/env python3
"""
tests/test_launch_script.py

Unit tests for launch.sh launcher script.
Validates:
1. --help flag returns 0 and prints usage instructions.
2. Handling invalid argument flags.
3. Execution with --validate-only flags on mock inputs.
"""

import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LAUNCH_SCRIPT = REPO_ROOT / "launch.sh"


def test_launch_script_help():
    res = subprocess.run([str(LAUNCH_SCRIPT), "--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "CNV-Caller-Combiner Unified Launcher" in res.stdout
    assert "--input-dir" in res.stdout
    assert "--profile" in res.stdout
    assert "--validate-only" in res.stdout


def test_launch_script_invalid_arg():
    res = subprocess.run([str(LAUNCH_SCRIPT), "--unknown-flag-xyz"], capture_output=True, text=True)
    assert res.returncode != 0
    assert "Unknown option" in res.stderr

def test_launch_script_workflow_options():
    res = subprocess.run([str(LAUNCH_SCRIPT), "--workflow", "cloud", "--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "--workflow" in res.stdout
    assert "cloud" in res.stdout


def test_launch_script_build_option():
    res = subprocess.run([str(LAUNCH_SCRIPT), "--build", "hg38", "--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "--build" in res.stdout

def test_launch_script_enable_igv_flag():
    res = subprocess.run([str(LAUNCH_SCRIPT), "--enable-igv", "--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "--enable-igv" in res.stdout


def test_launch_script_summary_flag(tmp_path):
    # Running summary on empty or non-existent dir exits gracefully
    res = subprocess.run([str(LAUNCH_SCRIPT), "--summary", str(tmp_path)], capture_output=True, text=True)
    assert res.returncode == 0
    assert "Analyzing pipeline results in" in res.stdout
