#!/usr/bin/env python3
"""Validate CI configuration files used by GitHub Actions and CircleCI."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "tests.yml"
CODEQL_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "codeql.yml"
CIRCLECI_CONFIG = REPO_ROOT / ".circleci" / "config.yml"


def test_tests_workflow_runs_repo_test_command():
    text = TESTS_WORKFLOW.read_text(encoding="utf-8")
    assert "python -m pytest -q" in text


def test_tests_workflow_checks_all_bash_shebang_scripts():
    text = TESTS_WORKFLOW.read_text(encoding="utf-8")
    assert "for f in bin/*; do" in text
    assert "head -n 1 \"$f\" | grep -Eq '^#!/(usr/bin/env )?bash'" in text


def test_tests_workflow_parses_r_scripts():
    text = TESTS_WORKFLOW.read_text(encoding="utf-8")
    assert "Parse R scripts" in text
    assert "Rscript --vanilla -e \"parse(file='$f')\"" in text


def test_codeql_workflow_includes_expected_languages():
    text = CODEQL_WORKFLOW.read_text(encoding="utf-8")
    assert "- language: actions" in text
    assert "- language: python" in text
    assert "- language: c-cpp" in text
    assert text.count("build-mode: none") >= 3


def test_tests_workflow_checks_nextflow_syntax():
    text = TESTS_WORKFLOW.read_text(encoding="utf-8")
    assert "nextflow" in text.lower(), (
        "tests.yml must install and run Nextflow to validate .nf file syntax"
    )
    assert "main.nf" in text, (
        "tests.yml must run 'nextflow run main.nf' (e.g. with --help) "
        "to catch Nextflow/Groovy syntax errors in workflow files"
    )


def test_circleci_config_runs_python_tests():
    text = CIRCLECI_CONFIG.read_text(encoding="utf-8")
    assert "test:" in text
    assert "cimg/python:3.10" in text
    assert "python -m pytest -q" in text
