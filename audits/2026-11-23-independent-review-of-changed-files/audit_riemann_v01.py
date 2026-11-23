#!/usr/bin/env python3
"""Non-destructive, independent v0.1 repository audit driver.

Usage:
    python audit_riemann_v01.py /path/to/riemann-repository
    python audit_riemann_v01.py . --output /tmp/riemann-audit.json

Run in an isolated environment with the repository's dependencies installed.
This script does not install packages, intentionally edit source files, or certify a theorem.
Repository tests may create caches or experiment artifacts; inspect git status afterward.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import tempfile


class Audit:
    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.findings: list[dict[str, object]] = []

    def record(self, level: str, name: str, detail: str) -> None:
        self.findings.append({"level": level, "check": name, "detail": detail})
        print(f"[{level:4}] {name}: {detail}")

    def file(self, name: str, *, essential: bool = True) -> Path:
        path = self.root / name
        self.record(
            "PASS" if path.is_file() else ("FAIL" if essential else "WARN"),
            f"file:{name}",
            "present" if path.is_file() else "missing",
        )
        return path

    def command(self, name: str, args: list[str], timeout: int = 180) -> None:
        start = time.monotonic()
        try:
            result = subprocess.run(
                args, cwd=self.root, capture_output=True, text=True,
                timeout=timeout, check=False,
                env={**os.environ, "PYTHONUTF8": "1", "PYTHONDONTWRITEBYTECODE": "1"},
            )
            output = (result.stdout + "\n" + result.stderr).strip()
            elapsed = time.monotonic() - start
            self.record(
                "PASS" if result.returncode == 0 else "FAIL",
                name,
                f"exit={result.returncode}; seconds={elapsed:.1f}; "
                f"output_tail={output[-2400:]!r}",
            )
        except subprocess.TimeoutExpired:
            self.record("FAIL", name, f"timed out after {timeout}s")
        except OSError as exc:
            self.record("FAIL", name, f"could not start: {exc}")

    def static_checks(self) -> None:
        essentials = [
            "pyproject.toml", ".github/workflows/ci.yml", "AGENTS.md",
            "CLAUDE.md", "requirements-lock.txt", ".python-version",
            "tools/verify_certificates.py", "tools/check_registry.py",
            "tools/check_english_only.py", "src/riemann_research/independent_check.py",
            "tests/test_historical_failures.py", "tests/test_independent_check.py",
            "tests/ORACLES.md", "docs/derivations/what-a-certificate-establishes.md",
        ]
        for item in essentials:
            self.file(item)
        for item in [".github/workflows/research.yml", "LICENSE", "CITATION.cff",
                     "docs/repository-settings.md", "docs/research-plan.md"]:
            self.file(item, essential=False)

        workflow = self.root / ".github/workflows/ci.yml"
        if workflow.is_file():
            data = workflow.read_text(encoding="utf-8")
            has_push = bool(re.search(r"(?m)^\s*push\s*:", data))
            has_pr = bool(re.search(r"(?m)^\s*pull_request\s*:", data))
            self.record("PASS" if has_push and has_pr else "FAIL", "ci:triggers",
                        f"push={has_push}, pull_request={has_pr}")
            self.record("PASS" if re.search(r"(?m)^\s*contents\s*:\s*read\b", data)
                        else "FAIL", "ci:read_only", "contents: read required")
            self.record("PASS" if "persist-credentials: false" in data else "WARN",
                        "ci:checkout_credentials", "checkout should not persist credentials")
            uses = re.findall(r"(?m)^\s*-\s*uses\s*:\s*([^\s#]+)", data)
            bad = [u for u in uses if not re.search(r"@[0-9a-f]{40}$", u)]
            self.record("PASS" if uses and not bad else "FAIL", "ci:action_pins",
                        f"uses={len(uses)}, non_SHA_pins={bad}")
            for marker in ["python -m pytest", "ruff check", "verify_certificates.py",
                           "check_registry.py --check", "check_english_only.py"]:
                self.record("PASS" if marker in data else "FAIL",
                            "ci:step:" + marker, "found" if marker in data else "missing")
            if re.search(r"(?m)^\s*branches\s*:\s*\[main\]", data):
                self.record("WARN", "ci:push_scope",
                            "push runs on main only; PR runs on opened/updated PRs")

        pin = self.root / ".python-version"
        if pin.is_file():
            expected = pin.read_text().strip()
            active = f"{sys.version_info.major}.{sys.version_info.minor}"
            self.record("PASS" if expected.startswith(active) else "WARN", "python:version",
                        f"repo={expected}, active={active}; run audit in intended environment")
        lock = self.root / "requirements-lock.txt"
        if lock.is_file():
            requirements = [s.strip() for s in lock.read_text().splitlines()
                            if s.strip() and not s.lstrip().startswith(("#", "-"))]
            unpinned = [s for s in requirements if "==" not in s and " @ " not in s]
            self.record("PASS" if not unpinned else "WARN", "dependencies:versions",
                        f"requirements={len(requirements)}, not_exactly_pinned={unpinned[:8]}")
            self.record("WARN", "dependencies:hashes",
                        "check whether lock has --hash entries; versions alone do not fix wheels/OS")

        cert_dir = self.root / "certificates"
        if cert_dir.exists():
            certs = list(cert_dir.rglob("*.json"))
            self.record("PASS" if len(certs) == 7 else "WARN", "certificates:count",
                        f"found={len(certs)} JSON files; author reported seven")
        else:
            self.record("WARN", "certificates:location", "not found at certificates/; inspect tree")

        independent = self.root / "src" / "riemann_research" / "independent_check.py"
        if independent.is_file():
            data = independent.read_text(encoding="utf-8")
            if "simpson" in data.lower():
                self.record("WARN", "math:simpson_exactness",
                            "Simpson's rule is exact for piecewise polynomial integrands ONLY "
                            "when its application respects every breakpoint and degree bounds. "
                            "Review the actual implementation and assumptions.")
            if re.search(r"\b(float|np\.float|numpy)\b", data):
                self.record("WARN", "math:exact_verifier",
                            "floating-point related token found; review claims of exactness")

        guide = self.root / "docs/derivations/what-a-certificate-establishes.md"
        if guide.is_file():
            s = guide.read_text(encoding="utf-8").lower()
            self.record("PASS" if ("rh" in s and "certificate" in s)
                        else "WARN", "math:proof_boundary",
                        "confirm certificate -> zeta theorem is NOT asserted as proved")

    def independent_reference_math(self) -> None:
        # These are intentionally separate from repository imports.
        # For the normalized centered box h=1_[-1/2,1/2], R=(1-|a|)_+.
        # Integral of a*(1-a) from 0 to 1 is 1/6, yielding C=4/3.
        c = Fraction(1) + 2 * (Fraction(1, 2) - Fraction(1, 3))
        self.record("PASS" if c == Fraction(4, 3) else "FAIL", "oracle:fejer",
                    f"exact baseline C={c}; arithmetic check, not repo implementation")
        # Tail example with widths 1/2 and 3/2, q_2=-1/2.
        # a=5/4 is beyond all cross-overlap supports, leaving q_2^2 self-overlap.
        tail = Fraction(1, 4) * (Fraction(3, 2) - Fraction(5, 4)) / Fraction(3, 2)**2
        self.record("PASS" if tail == Fraction(1, 36) else "FAIL", "oracle:positive_tail",
                    f"R(5/4)={tail} > 0 for a signed centered-box mixture")
        self.record("WARN", "math:unproven_lemma",
                    "The reported integral_1^2 F <= 8/3 is proof-pending. "
                    "Passing certificate checks cannot establish this analytic inequality.")

    def execute(self, *, timeout: int = 180) -> None:
        py = sys.executable
        self.command("lint:ruff", [py, "-m", "ruff", "check", "--no-cache", "."], timeout)
        self.command("check:english", [py, "tools/check_english_only.py"], timeout)
        self.command("check:certificates", [py, "tools/verify_certificates.py"], timeout)
        self.command("check:registry", [py, "tools/check_registry.py", "--check"], timeout)
        self.command("tests:historical-and-independent",
                     [py, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/test_historical_failures.py",
                      "tests/test_independent_check.py"], timeout)
        self.command("tests:all", [py, "-m", "pytest", "-q", "-p", "no:cacheprovider"], max(timeout, 300))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", nargs="?", default=".", help="unpacked repository path")
    parser.add_argument("--output", default=None, help="JSON report path (prefer outside repository)")
    parser.add_argument("--static-only", action="store_true", help="do not execute repo code")
    parser.add_argument("--timeout", type=int, default=180, help="timeout per command in seconds")
    args = parser.parse_args()
    root = Path(args.repo).expanduser().resolve()
    if not (root / "pyproject.toml").is_file():
        print(f"ERROR: {root} is not the v0.1 repository root (missing pyproject.toml)")
        return 2

    audit = Audit(root)
    print(f"Independent audit of {root}")
    audit.static_checks()
    audit.independent_reference_math()
    if not args.static_only:
        audit.execute(timeout=args.timeout)

    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                         capture_output=True, text=True, check=False)
    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "repo": str(root),
        "git_head": rev.stdout.strip() if rev.returncode == 0 else None,
        "python": sys.version,
        "platform": platform.platform(),
        "static_only": args.static_only,
        "findings": audit.findings,
        "summary": {k: sum(x["level"] == k for x in audit.findings)
                    for k in ("PASS", "WARN", "FAIL")},
        "scope": "Software/CI and reference-arithmetic audit; NOT a mathematical proof",
    }
    output = Path(args.output or (Path(tempfile.gettempdir()) / "riemann_audit_report.json")).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nSummary: {report['summary']}\nReport: {output}")
    return 1 if report["summary"]["FAIL"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
