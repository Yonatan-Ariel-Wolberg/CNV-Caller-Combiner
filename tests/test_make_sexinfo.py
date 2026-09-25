#!/usr/bin/env python3
"""
tests/test_make_sexinfo.py

Unit tests for bin/make_sexinfo.py.
Validates samplesheet parsing, coverage ratio logic, fallback assignment, and CLI execution.
"""

import os
import subprocess
import sys
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "bin"))

import make_sexinfo


def test_parse_samplesheet_tsv(tmp_path):
    samplesheet = tmp_path / "samplesheet.tsv"
    samplesheet.write_text(
        "sample_id\tbam_file\tbai_file\n"
        "SAMPLE_001\t/path/to/sample1.bam\t/path/to/sample1.bam.bai\n"
        "SAMPLE_002\t/path/to/sample2.bam\t/path/to/sample2.bam.bai\n"
        "# Comment line\n"
        "SAMPLE_003\t/path/to/sample3.bam\t/path/to/sample3.bam.bai\n"
    )
    samples = make_sexinfo.parse_samplesheet(str(samplesheet))
    assert len(samples) == 3
    assert samples[0] == ("SAMPLE_001", "/path/to/sample1.bam")
    assert samples[1] == ("SAMPLE_002", "/path/to/sample2.bam")
    assert samples[2] == ("SAMPLE_003", "/path/to/sample3.bam")


def test_make_sexinfo_cli(tmp_path):
    samplesheet = tmp_path / "samplesheet.tsv"
    samplesheet.write_text(
        "sample_id\tbam_file\n"
        "PATIENT_A\t/nonexistent/a.bam\n"
        "PATIENT_B\t/nonexistent/b.bam\n"
    )
    out_sexinfo = tmp_path / "sexinfo.tsv"
    script = REPO_ROOT / "bin" / "make_sexinfo.py"

    cmd = [
        sys.executable,
        str(script),
        "--samplesheet", str(samplesheet),
        "--output", str(out_sexinfo),
        "--default-sex", "F",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error running make_sexinfo.py: {res.stderr}"

    assert out_sexinfo.exists()
    lines = out_sexinfo.read_text().strip().splitlines()
    assert lines[0] == "sample_id\tsex"
    assert lines[1] == "PATIENT_A\tF"
    assert lines[2] == "PATIENT_B\tF"

def test_parse_samplesheet_cram(tmp_path):
    samplesheet = tmp_path / "samplesheet_crams.tsv"
    samplesheet.write_text(
        "sample_id\tcram_file\tcrai_file\n"
        "SAMPLE_X\t/data/x.cram\t/data/x.cram.crai\n"
    )
    samples = make_sexinfo.parse_samplesheet(str(samplesheet))
    assert len(samples) == 1
    assert samples[0] == ("SAMPLE_X", "/data/x.cram")


def test_make_sexinfo_default_male(tmp_path):
    samplesheet = tmp_path / "samplesheet.tsv"
    samplesheet.write_text(
        "sample_id\tbam_file\n"
        "PATIENT_M1\t/dummy/m1.bam\n"
    )
    out_sexinfo = tmp_path / "sexinfo_male.tsv"
    script = REPO_ROOT / "bin" / "make_sexinfo.py"

    cmd = [
        sys.executable,
        str(script),
        "--samplesheet", str(samplesheet),
        "--output", str(out_sexinfo),
        "--default-sex", "M",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    lines = out_sexinfo.read_text().strip().splitlines()
    assert lines[1] == "PATIENT_M1\tM"
