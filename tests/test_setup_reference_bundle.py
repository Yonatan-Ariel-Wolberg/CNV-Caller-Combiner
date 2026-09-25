#!/usr/bin/env python3
"""
tests/test_setup_reference_bundle.py

Unit tests for bin/setup_reference_bundle.sh.
Validates help command and staging with --skip-fasta into a temporary directory.
"""

import json
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "bin" / "setup_reference_bundle.sh"


def test_setup_reference_bundle_help():
    res = subprocess.run([str(SCRIPT), "--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "CNV-CALLER-COMBINER: REFERENCE BUNDLE SETUP" in res.stdout or "setup_reference_bundle.sh" in res.stdout
    assert "--skip-fasta" in res.stdout


def test_setup_reference_bundle_skip_fasta(tmp_path):
    dest_dir = tmp_path / "test_refs"
    res = subprocess.run([str(SCRIPT), "--dest", str(dest_dir), "--skip-fasta"], capture_output=True, text=True)
    assert res.returncode == 0

    json_path = dest_dir / "reference_paths.json"
    assert json_path.exists()
    data = json.loads(json_path.read_text())
    assert "ref" in data
    assert "xhmm_conf" in data
    assert (dest_dir / "xhmm" / "params.txt").exists()
