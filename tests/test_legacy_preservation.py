"""The pre-audit code is preserved byte for byte and still reproduces its archived outputs.

No translation was needed (the legacy tree contains no non-English text), so these tests also
serve as the 'translation introduced no behavioural change' regression check required by the
English-only policy.
"""
import csv
import hashlib
import json

import numpy as np
import pytest


def test_legacy_files_match_the_recorded_sha256_manifest(root):
    base = root / "legacy" / "pair_correlation_sensitivity_v0"
    lines = (root / "legacy" / "MANIFEST.sha256").read_text().split("\n")
    entries = [line.split(None, 1) for line in lines if line.strip()]
    assert len(entries) == 5
    for digest, rel in entries:
        assert hashlib.sha256((base / rel.strip()).read_bytes()).hexdigest() == digest, rel


def test_legacy_positive_scan_reproduces_archived_summary(root):
    import pair_sensitivity as legacy_pos
    archived = json.loads((root / "audits" / "2026-11-23-pair-correlation-audit" /
                           "legacy_run_outputs" / "summary.json").read_text())
    for A in [1.1, 1.5, 2.0]:
        stages = legacy_pos.adaptive_scan(A, refinement_rounds=2, top_k=3, basis_count=12)
        for s in stages:
            ref = next(r for r in archived if r["A"] == A and r["stage"] == s["stage"])
            assert s["C"] == pytest.approx(ref["C_candidate_conditional"], abs=1e-9)
            assert [int(i) for i in s["top_bins"]] == ref["top_bins"]


def test_legacy_signed_scan_reproduces_archived_values(root):
    import signed_sensitivity as legacy_signed
    path = root / "audits" / "2026-11-23-pair-correlation-audit" / "legacy_run_outputs" / "signed_scan.csv"
    rows = list(csv.DictReader(path.open()))
    edges = np.linspace(1, 2, 11)
    for row in rows:
        A, s = float(row["A"]), float(row["cap_scale"])
        if (A, s) not in [(1.5, 4.0), (2.0, 4.0)]:
            continue
        pr = legacy_signed.make_program(A, 12, edges, masses=s * np.diff(edges))
        signed = legacy_signed.optimize_signed(pr, mode="caps", restarts=3, seed=42)
        assert signed["P"] == pytest.approx(float(row["P_signed"]), abs=1e-9)
