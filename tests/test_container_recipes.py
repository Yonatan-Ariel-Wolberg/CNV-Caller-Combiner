#!/usr/bin/env python3
"""Validate container recipe files for ICAv2 and train workflows."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
ICAV2_DOCKERFILE = REPO_ROOT / "bin" / "Dockerfile.icav2"
TRAIN_DOCKERFILE = REPO_ROOT / "bin" / "Dockerfile.train"
README = REPO_ROOT / "README.md"


def test_icav2_dockerfile_exists_and_has_expected_content():
    text = ICAV2_DOCKERFILE.read_text(encoding="utf-8")
    assert "FROM ubuntu:24.04" in text
    assert "software.version=\"2.43.0\"" in text
    assert "ica-linux-amd64.zip" in text
    assert "/usr/local/bin/icav2" in text


def test_train_dockerfile_exists_and_has_expected_content():
    text = TRAIN_DOCKERFILE.read_text(encoding="utf-8")
    assert "FROM quay.io/condaforge/mambaforge:24.9.2-0" in text
    assert "xgboost=2.1.4" in text
    assert "scikit-learn=1.6.1" in text
    assert "imbalanced-learn=0.13.0" in text
    assert "pandas=2.2.3" in text
    assert "numpy=2.2.3" in text


def test_readme_documents_docker_build_commands_for_new_recipes():
    text = README.read_text(encoding="utf-8")
    assert "docker build -f bin/Dockerfile.icav2" in text
    assert "docker build -f bin/Dockerfile.train" in text
