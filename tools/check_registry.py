"""Validate results/registry.yaml against results/assumptions.yaml and the certificates.

Usage:
    python tools/check_registry.py            # validate and print the generated table
    python tools/check_registry.py --write    # also regenerate docs/results-table.md
    python tools/check_registry.py --check    # fail if docs/results-table.md is out of date

Exit status 0 means the registry is internally consistent and every certificate passes BOTH exact
verifiers. It never assigns or upgrades a mathematical status, and it does not mean any result is proven.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from riemann_research import certify, independent_check  # noqa: E402
from riemann_research.status import (  # noqa: E402
    DISCLAIMER, Assumption, Computation, Mathematical, Optimality, conditional_on, consistency_errors)

KINDS = {"result", "negative-result", "known-defect"}
FORBIDDEN_FIELDS = {"label", "status", "verified", "proven", "argument_status"}
PROVENANCE_FIELDS = {"commit", "generator", "seed", "precision", "solver", "verifiers"}
TABLE = ROOT / "docs" / "results-table.md"
COMMON = ("montgomery-theorem", "F-nonnegative", "simple-zero-inequality")


class RegistryError(Exception):
    pass


def load():
    reg = yaml.safe_load((ROOT / "results" / "registry.yaml").read_text())
    ass = yaml.safe_load((ROOT / "results" / "assumptions.yaml").read_text())
    return reg, ass


def _cap_scope_errors(eid, cert, cited, ass):
    """Every bin that can carry mass must be covered by the scope of a cited cap assumption."""
    out = cert["outside"]
    if out["mode"] != "caps":
        return []
    scoped = [ass[a]["scope"] for a in cited if "scope" in ass[a]]
    if not scoped:
        return [f"{eid}: caps-mode certificate cites no cap assumption with a scope"]
    S = max(Fraction(w) for w in cert["widths"])
    edges = [Fraction(e) for e in out["edges"]]
    errs = []
    for (lo, hi), cap in zip(zip(edges[:-1], edges[1:], strict=True), map(Fraction, out["caps"]), strict=True):
        if lo >= S:
            continue
        ok = False
        for sc in scoped:
            if "max_bin_length" in sc and (hi - lo > Fraction(sc["max_bin_length"]) or cap < Fraction(sc["min_cap"])):
                continue
            if "cap_per_unit_width" in sc and (cap < Fraction(sc["cap_per_unit_width"]) * (hi - lo)
                                               or hi > Fraction(sc["max_edge"])):
                continue
            ok = True
            break
        if not ok:
            errs.append(f"{eid}: bin [{float(lo)}, {float(hi)}] with cap {float(cap)} is outside the scope of "
                        "every cited cap assumption")
    return errs


@functools.lru_cache(maxsize=None)
def _verify_both(content_sha256: str, text: str):
    """Run both exact verifiers once per distinct certificate content (cache keyed by content hash)."""
    cert = json.loads(text)
    return certify.verify(cert), independent_check.check(cert)


def _check_certificate(eid, e, ass, errors):
    path = ROOT / e.get("certificate", "")
    if not e.get("certificate") or not path.is_file():
        errors.append(f"{eid}: numerically-certified requires an existing certificate file")
        return
    text = path.read_text()
    cert = json.loads(text)
    missing = PROVENANCE_FIELDS - set(cert.get("provenance", {}))
    if missing:
        errors.append(f"{eid}: certificate provenance lacks {sorted(missing)}")
    try:
        a, b = _verify_both(hashlib.sha256(text.encode()).hexdigest(), text)
    except (certify.CertificateError, ValueError) as exc:
        errors.append(f"{eid}: certificate failed: {exc}")
        return
    if any(a[k] != b[k] for k in ("C", "P", "outside_term", "max_R_on_1_S")):
        errors.append(f"{eid}: the two independent verifiers disagree")
    if Fraction(str(e["value"])) > a["P"]:
        errors.append(f"{eid}: registry value {e['value']} exceeds the certified P {float(a['P'])}")
    if abs(float(a["P"]) - float(e["value"])) > 1e-9:
        errors.append(f"{eid}: registry value {e['value']} differs from certified P {float(a['P']):.12f}")
    errors.extend(_cap_scope_errors(eid, cert, e["assumptions"], ass))


def validate(reg: dict, ass: dict) -> list[dict]:
    errors = []
    if reg.get("schema") != "riemann-research/registry/v1":
        errors.append("unknown registry schema")
    for aid, a in ass.items():
        try:
            Assumption(a["status"])
        except (KeyError, ValueError):
            errors.append(f"assumption {aid}: invalid status {a.get('status')!r}")
    record = reg["literature"]["chirre-goncalves-de-laat-2020"]["value"]
    ids = [e.get("id") for e in reg["results"]]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        errors.append(f"duplicate ids: {sorted(dup)}")
    rows = []
    for e in reg["results"]:
        eid = e.get("id", "?")
        kind = e.get("kind")
        if kind not in KINDS:
            errors.append(f"{eid}: unknown kind {kind!r}")
            continue
        bad = FORBIDDEN_FIELDS & set(e)
        if bad:
            errors.append(f"{eid}: hand-written or obsolete status fields are not allowed: {sorted(bad)}")
        if kind == "known-defect":
            for t in e.get("tests", []):
                if not (ROOT / t).is_file():
                    errors.append(f"{eid}: test file {t} missing")
            continue
        try:
            comp = Computation(e["computation"])
            math = Mathematical(e["mathematical_status"])
        except (KeyError, ValueError) as exc:
            errors.append(f"{eid}: missing or invalid status: {exc}")
            continue
        arguments = e.get("argument", [])
        for ref in arguments:
            if not (ROOT / ref).is_file():
                errors.append(f"{eid}: argument {ref} missing")
        errors.extend(f"{eid}: {msg}" for msg in consistency_errors(
            comp, math, e.get("reviews", []), e.get("open_issues", []), bool(arguments),
            numeric_claim=(kind == "result")))
        if kind == "negative-result":
            for ref in e.get("evidence", []) + e.get("tests", []):
                if ref not in ids and not (ROOT / ref).exists():
                    errors.append(f"{eid}: evidence {ref!r} is neither a registry id nor a file")
            continue
        try:
            Optimality(e["optimality"])
        except (KeyError, ValueError) as exc:
            errors.append(f"{eid}: {exc}")
        missing = [a for a in e.get("assumptions", []) if a not in ass]
        if missing:
            errors.append(f"{eid}: unknown assumptions {missing}")
            continue
        if "RH" not in e["assumptions"]:
            errors.append(f"{eid}: every result in this repository must list RH")
        if comp is Computation.NUMERICALLY_CERTIFIED:
            _check_certificate(eid, e, ass, errors)
        elif e.get("certificate"):
            errors.append(f"{eid}: a certificate is attached but computation is {comp.value}")
        statuses = [Assumption(ass[a]["status"]) for a in e["assumptions"]]
        # The pair-correlation benchmark is not an upper bound. A value above it is allowed, but only
        # with an explicit, human-written flag so that it is audited before anyone reports it.
        above = (conditional_on(statuses) != "RH + UNPROVEN hypotheses" and float(e["value"]) > record)
        flag = e.get("record_comparison", {})
        if above and not (flag.get("flagged_for_audit") is True and flag.get("note")):
            errors.append(f"{eid}: value {e['value']} exceeds the pair-correlation benchmark {record}; add "
                          "'record_comparison: {flagged_for_audit: true, note: ...}' after a deliberate check")
        marker = " [EXCEEDS PAIR-CORRELATION BENCHMARK - flagged for audit]" if above else ""
        rows.append({"id": eid, "value": e["value"], "computation": comp.value, "mathematical": math.value,
                     "optimality": e["optimality"], "conditional_on": conditional_on(statuses),
                     "assumptions": ", ".join(a for a in e["assumptions"] if a not in COMMON and a != "RH"),
                     "title": e["title"] + marker})
    if errors:
        raise RegistryError("\n".join(errors))
    return rows


def render(rows, reg) -> str:
    lines = ["# Results table", "",
             "<!-- Generated by tools/check_registry.py --write. Do not edit by hand. -->", "",
             f"> {DISCLAIMER}", "",
             "P is a lower bound for the proportion of simple zeros. Every row assumes RH and Montgomery's theorem.",
             "No row is `proven` or `independently-reviewed`: those statuses require a recorded human review.",
             "The pair-correlation benchmark under RH is "
             f"{reg['literature']['chirre-goncalves-de-laat-2020']['value']} (Chirre-Goncalves-de Laat 2020); "
             "it is not the best known proportion of simple zeros by all methods.", "",
             "| Mathematical status | Computation | Conditional on | P | Extra assumptions | Result |",
             "|---|---|---|---|---|---|"]
    rank = {"RH": 0, "RH + unreviewed lemma": 1, "RH + UNPROVEN hypotheses": 2}
    for r in sorted(rows, key=lambda r: (rank[r["conditional_on"]], r["mathematical"] != "proof-pending",
                                         -float(r["value"]))):
        lines.append(f"| {r['mathematical']} | {r['computation']} | {r['conditional_on']} | {r['value']} | "
                     f"{r['assumptions'] or '-'} | `{r['id']}`: {r['title']} |")
    lines += ["", "Negative results and known defects are listed in `results/registry.yaml`.", ""]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    reg, ass = load()
    try:
        rows = validate(reg, ass)
    except RegistryError as exc:
        print("REGISTRY CHECK FAILED\n" + str(exc))
        return 1
    text = render(rows, reg)
    if args.write:
        TABLE.write_text(text)
    if args.check and (not TABLE.is_file() or TABLE.read_text() != text):
        print("docs/results-table.md is out of date; run: python tools/check_registry.py --write")
        return 1
    print(text)
    print(f"Registry consistent: {len(rows)} results; every certificate passed both exact verifiers.")
    print("Note: " + DISCLAIMER)
    return 0


if __name__ == "__main__":
    sys.exit(main())
