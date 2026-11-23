"""Experiment metadata recorded with every generated certificate or output.

Fields: repository commit (and whether the tree was dirty), generator script, inputs, seed,
precision (rounding denominator), solver and its version, verifier versions, library versions,
and a UTC timestamp. Metadata is informational; verifiers never use it to decide validity.
"""
from __future__ import annotations

import datetime as _dt
import platform
import subprocess
from importlib import metadata
from pathlib import Path

from . import __version__

VERIFIERS = {"certify": "riemann-research/box-certificate/v1", "independent_check": "v1",
             "package": __version__}


def _version(dist: str) -> str:
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return "not installed"


def _git(root: Path, *args: str) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True,
                              check=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def collect(generator: str, *, inputs: dict | None = None, seed=None, precision=None,
            solver: str | None = None, solver_package: str | None = None, root: Path | None = None) -> dict:
    root = Path(root) if root else Path(__file__).resolve().parents[2]
    commit = _git(root, "rev-parse", "HEAD")
    dirty = _git(root, "status", "--porcelain") if commit else None
    solver_version = _version(solver_package) if solver_package else None
    return {
        "commit": commit or "unavailable (not a git checkout)",
        "dirty": bool(dirty) if commit else None,
        "generator": generator,
        "inputs": inputs or {},
        "seed": seed,
        "precision": precision,
        "solver": solver,
        "solver_version": solver_version,
        "verifiers": VERIFIERS,
        "python": platform.python_version(),
        "libraries": {d: _version(d) for d in ("numpy", "scipy", "cvxpy")},
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
