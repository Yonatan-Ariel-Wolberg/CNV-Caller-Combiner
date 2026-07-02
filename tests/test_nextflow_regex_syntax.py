#!/usr/bin/env python3
"""Guard against regex syntax and Groovy patterns that break Nextflow parsing."""

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

# file(params.xxx ?: 'NO_FILE') — the Elvis operator and fallback value fit
# comfortably within 80 characters of the opening 'file(' token.
_ELVIS_SEARCH_WINDOW = 80


def test_nf_files_do_not_use_trailing_regex_i_flag():
    def has_trailing_i_flag(text: str) -> bool:
        idx = text.find("/i")
        while idx != -1:
            next_pos = idx + 2
            if next_pos >= len(text):
                return True
            next_char = text[next_pos]
            if not (next_char.isalnum() or next_char == "_"):
                return True
            idx = text.find("/i", idx + 2)
        return False

    for nf_path in REPO_ROOT.rglob("*.nf"):
        text = nf_path.read_text(encoding="utf-8")
        assert not has_trailing_i_flag(text), (
            f"{nf_path} contains /.../i regex syntax; use inline (?i) instead"
        )


def test_nf_files_do_not_call_file_with_empty_param_default():
    """file(params.xxx) at module level must use ?: 'NO_FILE' if param defaults to ''.

    Nextflow raises "Argument of `file()` function cannot be empty" when
    file() receives an empty string.  Any module that sets ``params.xxx = ''``
    and then calls ``file(params.xxx)`` without the Groovy Elvis fallback
    (``params.xxx ?: 'NO_FILE'``) will fail the ``nextflow run --help`` syntax
    check in CI.
    """
    empty_default_re = re.compile(r"^params\.(\w+)\s*=\s*''", re.MULTILINE)
    file_call_re = re.compile(r"\bfile\(\s*params\.(\w+)")

    for nf_path in REPO_ROOT.rglob("*.nf"):
        text = nf_path.read_text(encoding="utf-8")
        empty_params = {m.group(1) for m in empty_default_re.finditer(text)}
        if not empty_params:
            continue
        for m in file_call_re.finditer(text):
            param_name = m.group(1)
            if param_name not in empty_params:
                continue
            # The call must use the ?: (Groovy Elvis) fallback operator so that
            # an empty default never reaches file().
            call_chunk = text[m.start(): m.start() + _ELVIS_SEARCH_WINDOW]
            assert "?:" in call_chunk, (
                f"{nf_path}: file(params.{param_name}) is called but "
                f"params.{param_name} defaults to '' — "
                f"use file(params.{param_name} ?: 'NO_FILE') to avoid "
                "'Argument of `file()` function cannot be empty' in CI"
            )
