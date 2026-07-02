#!/usr/bin/env python3
"""Docs tests for module-specific documentation split."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
README = REPO_ROOT / "README.md"
DOCS_DIR = REPO_ROOT / "docs" / "modules"


MODULE_DOCS = {
    "canoes": "canoes.md",
    "xhmm": "xhmm.md",
    "clamms": "clamms.md",
    "gcnv": "gatk-gcnv.md",
    "cnvkit": "cnvkit.md",
    "dragen": "dragen.md",
    "indelible": "indelible.md",
    "survivor": "survivor.md",
    "truvari": "truvari.md",
    "feature_extraction": "feature-extraction.md",
    "train": "train.md",
    "evaluate": "evaluate.md",
    "normalise": "normalise.md",
}


def test_readme_links_to_each_module_doc():
    readme_text = README.read_text(encoding="utf-8")
    for doc_file in MODULE_DOCS.values():
        rel_link = f"docs/modules/{doc_file}"
        assert rel_link in readme_text, f"README is missing link to {rel_link}"


def test_each_module_doc_contains_required_sections():
    for workflow, doc_file in MODULE_DOCS.items():
        text = (DOCS_DIR / doc_file).read_text(encoding="utf-8")
        assert "## Methodology" in text, f"{doc_file} missing Methodology section"
        assert "## Processes in this module" in text, f"{doc_file} missing process section"
        assert "## Required parameters" in text, f"{doc_file} missing required parameters section"
        assert f"`workflow`: `{workflow}`" in text, f"{doc_file} missing workflow parameter"
