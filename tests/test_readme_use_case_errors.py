#!/usr/bin/env python3
"""Ensure README lists currently unhandled use-case validation gaps."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
README = REPO_ROOT / "README.md"


def test_readme_documents_validation_gaps():
    text = README.read_text(encoding="utf-8")
    assert "### Known unhandled use-case errors (current gaps)" in text

    expected_items = [
        "`--workflow gcnv`",
        "`--workflow cnvkit`",
        "`--workflow feature_extraction`",
        "`--workflow train`",
        "Any workflow: some workflow-specific secondary assets may exist but still be",
        "unreadable or invalid in ways not fully pre-validated.",
    ]
    for item in expected_items:
        assert item in text, f"README missing unhandled use-case error bullet: {item}"
