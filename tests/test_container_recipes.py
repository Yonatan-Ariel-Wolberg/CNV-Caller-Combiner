#!/usr/bin/env python3
"""Validate container recipe files for ICAv2 and train workflows."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
ICAV2_DOCKERFILE = REPO_ROOT / "containers" / "Dockerfile.icav2"
TRAIN_DOCKERFILE = REPO_ROOT / "containers" / "Dockerfile.train"
SEENV_DOCKERFILE = REPO_ROOT / "containers" / "Dockerfile.seenv"
BIN_ICAV2_DOCKERFILE = REPO_ROOT / "bin" / "Dockerfile.icav2"
BIN_TRAIN_DOCKERFILE = REPO_ROOT / "bin" / "Dockerfile.train"
BIN_SEENV_DOCKERFILE = REPO_ROOT / "bin" / "Dockerfile.seenv"
README = REPO_ROOT / "README.md"


def test_icav2_dockerfile_exists_and_has_expected_content():
    text = ICAV2_DOCKERFILE.read_text(encoding="utf-8")
    assert "FROM ubuntu:24.04" in text
    assert "software.version=\"2.43.0\"" in text
    assert "ica-linux-amd64.zip" in text
    assert "/usr/local/bin/icav2" in text


def test_train_dockerfile_exists_and_has_expected_content():
    text = TRAIN_DOCKERFILE.read_text(encoding="utf-8")
    assert "FROM python:3.12-slim" in text or "FROM quay.io/condaforge/mambaforge" in text
    assert 'org.opencontainers.image.version="1.0.0"' in text
    assert "catboost" in text
    assert "shap" in text
    assert "scikit-learn" in text
    assert "imbalanced-learn" in text
    assert "pandas" in text
    assert "numpy" in text


def test_seenv_dockerfile_exists_and_has_expected_content():
    text = SEENV_DOCKERFILE.read_text(encoding="utf-8")
    assert "FROM quay.io/condaforge/mambaforge:24.9.2-0" in text
    assert "samtools=1.21" in text
    assert "bedtools=2.31.1" in text
    assert "mosdepth=0.3.14" in text
    assert "pysam=0.22.1" in text
    assert "matplotlib=3.8.4" in text
    assert "pandas=2.2.3" in text
    assert "numpy=1.26.4" in text
    assert "scipy=1.13.1" in text
    assert "MSBradshaw/SeeNV" in text
    assert "95566019d1ac880efe1b658623d74b29314376bc" in text
    assert "USER 10001:10001" in text
    assert "05638fb74bca9a1a0a888b70b084b039771b61362b997bab73487ee0635c1171" in text


def test_readme_documents_docker_build_commands_for_new_recipes():
    text = README.read_text(encoding="utf-8")
    assert "docker build -f bin/Dockerfile.icav2" in text or "containers/Dockerfile.icav2" in text
    assert "docker build -f bin/Dockerfile.train" in text or "containers/Dockerfile.train" in text
    assert "docker build -f bin/Dockerfile.seenv" in text or "containers/Dockerfile.seenv" in text


def test_bin_symlinks_exist_and_target_containers():
    assert BIN_ICAV2_DOCKERFILE.is_symlink()
    assert BIN_TRAIN_DOCKERFILE.is_symlink()
    assert BIN_SEENV_DOCKERFILE.is_symlink()
    assert BIN_ICAV2_DOCKERFILE.resolve() == ICAV2_DOCKERFILE.resolve()
    assert BIN_TRAIN_DOCKERFILE.resolve() == TRAIN_DOCKERFILE.resolve()
    assert BIN_SEENV_DOCKERFILE.resolve() == SEENV_DOCKERFILE.resolve()


