#!/usr/bin/env python3
"""Python entry point that reuses the log_message.py shipped with thisutils.

    from log_message import log_message          # this directory on sys.path
    log_message("message", message_type="info")

Resolution order:
  1. $LOG_MESSAGE_PY (explicit path to thisutils/log_message.py)
  2. R lookup: log_message.py under Rscript -e 'cat(system.file("scripts", package="thisutils"))'
  3. Plain fallback with the same call signature and ignored options

Progress and information messages only; stdout parsed by other scripts (for
example CASE_DONE lines) is left untouched.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

_cache: dict[str, object] = {}


def _thisutils_python_dir() -> Path | None:
    explicit = os.environ.get("LOG_MESSAGE_PY")
    if explicit and Path(explicit).is_file():
        return Path(explicit).resolve().parent
    # Current layout: thisutils/scripts/; earlier installs used python/.
    expr = ('d <- system.file("scripts", package = "thisutils"); '
            'if (!nzchar(d)) d <- system.file("python", package = "thisutils"); cat(d)')
    try:
        out = subprocess.run([os.environ.get("BRAINOMICS_RSCRIPT", "Rscript"), "--vanilla", "-e", expr],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    candidate = out.stdout.strip()
    if out.returncode == 0 and candidate and Path(candidate, "log_message.py").is_file():
        return Path(candidate)
    return None


def _load():
    if "fn" in _cache:
        return _cache["fn"]
    directory = _thisutils_python_dir()
    if directory is not None:
        # Load by path: importing log_message here would recurse into this module.
        import importlib.util  # noqa: PLC0415

        module_path = directory / "log_message.py"
        spec = importlib.util.spec_from_file_location("_thisutils_log_message", module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _cache["fn"] = module.log_message
        _cache["module"] = module
        return module.log_message

    def fallback(message, *args, **kwargs):
        options = {"message_type", "timestamp", "verbose", "cli_model", "indent",
                   "color", "bg_color", "multiline_indent"}
        extra = " ".join(str(a) for a in args if a not in options)
        print(message if not extra else f"{message} {extra}", file=sys.stderr)

    print("log_message: thisutils not found, falling back to plain output "
          "(set LOG_MESSAGE_PY)", file=sys.stderr)
    _cache["fn"] = fallback
    return fallback


def log_message(message, *args, **kwargs):  # noqa: D103
    return _load()(message, *args, **kwargs)


__all__ = ["log_message"]
