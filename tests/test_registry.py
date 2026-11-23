"""Registry consistency: three status dimensions, cap scopes, and honest comparisons."""
import json

import pytest
import yaml

import check_registry


@pytest.fixture(scope="module")
def loaded():
    reg, ass = check_registry.load()
    return reg, ass, check_registry.validate(reg, ass)


def copy(reg):
    return yaml.safe_load(yaml.safe_dump(reg))


def entry(reg, eid):
    return next(e for e in reg["results"] if e["id"] == eid)


def test_registry_validates_and_all_certificates_pass_both_verifiers(loaded):
    _, _, rows = loaded
    assert len(rows) >= 20


def test_generated_results_table_is_up_to_date(loaded, root):
    reg, _, rows = loaded
    assert (root / "docs" / "results-table.md").read_text() == check_registry.render(rows, reg)


def test_nothing_in_the_registry_is_reviewed_or_proven_yet(loaded):
    reg, _, _ = loaded
    for e in reg["results"]:
        if e["kind"] != "known-defect":
            assert e["mathematical_status"] in ("conjectural", "proof-pending")
            assert not e.get("reviews")


def test_results_with_unproven_hypotheses_are_flagged_in_the_table(loaded):
    reg, ass, rows = loaded
    by_id = {r["id"]: r for r in rows}
    for e in reg["results"]:
        if e["kind"] == "result" and any(ass[a]["status"] == "unproven-hypothesis" for a in e["assumptions"]):
            assert by_id[e["id"]]["conditional_on"] == "RH + UNPROVEN hypotheses"


def test_the_71_percent_figure_is_certified_arithmetic_but_conditional(loaded):
    _, _, rows = loaded
    r = next(r for r in rows if r["id"] == "pc-legacy-caps-A2-scale4")
    assert r["computation"] == "numerically-certified"
    assert r["mathematical"] == "proof-pending"
    assert r["conditional_on"] == "RH + UNPROVEN hypotheses"


def test_values_above_the_pair_correlation_benchmark_need_an_explicit_audit_flag(loaded):
    """The benchmark is not an upper bound: exceeding it is allowed but must be flagged by a person."""
    reg, ass, _ = loaded
    bad = copy(reg)
    bad["literature"]["chirre-goncalves-de-laat-2020"]["value"] = 0.675   # pretend the benchmark were lower
    with pytest.raises(check_registry.RegistryError, match="exceeds the pair-correlation benchmark"):
        check_registry.validate(bad, ass)
    for e in bad["results"]:
        if e["id"] in ("pc-provencap-A2", "pc-baseline-A1") or (
                e["kind"] == "result" and e["id"].startswith(("pc-legacy-signonly", "pc-provencap"))):
            e["record_comparison"] = {"flagged_for_audit": True, "note": "deliberately checked"}
    rows = check_registry.validate(bad, ass)
    flagged = [r for r in rows if "EXCEEDS PAIR-CORRELATION BENCHMARK" in r["title"]]
    assert {r["id"] for r in flagged} == {"pc-provencap-A2"}
    unproven = [r for r in rows if r["conditional_on"] == "RH + UNPROVEN hypotheses"]
    assert all("EXCEEDS" not in r["title"] for r in unproven)   # conditional values are not compared


@pytest.mark.parametrize("status", ["independently-reviewed", "proven"])
def test_upgrading_without_a_human_review_is_refused(loaded, status):
    reg, ass, _ = loaded
    bad = copy(reg)
    e = entry(bad, "pc-baseline-A1")
    e["mathematical_status"] = status
    e["reviews"] = [{"kind": "ai-audit", "reviewer": "Claude", "date": "2026-11-23", "scope": "all"}]
    with pytest.raises(check_registry.RegistryError, match="human review"):
        check_registry.validate(bad, ass)


def test_numerical_results_cannot_be_proof_pending(loaded):
    reg, ass, _ = loaded
    bad = copy(reg)
    entry(bad, "pc-legacy-caps-A1.1-scale1")["mathematical_status"] = "proof-pending"
    with pytest.raises(check_registry.RegistryError, match="numerically certified"):
        check_registry.validate(bad, ass)


def test_hand_written_labels_are_refused(loaded):
    reg, ass, _ = loaded
    bad = copy(reg)
    bad["results"][0]["label"] = "Verified"
    with pytest.raises(check_registry.RegistryError, match="hand-written"):
        check_registry.validate(bad, ass)


def test_certified_status_without_a_certificate_is_refused(loaded):
    reg, ass, _ = loaded
    bad = copy(reg)
    entry(bad, "pc-legacy-caps-A1.1-scale1")["computation"] = "numerically-certified"
    with pytest.raises(check_registry.RegistryError, match="certificate"):
        check_registry.validate(bad, ass)


def test_registry_value_above_the_certified_value_is_refused(loaded):
    reg, ass, _ = loaded
    bad = copy(reg)
    entry(bad, "pc-provencap-A2")["value"] = 0.6792
    with pytest.raises(check_registry.RegistryError, match="exceeds"):
        check_registry.validate(bad, ass)


def test_cap_used_outside_its_scope_is_refused(loaded, root, tmp_path, monkeypatch):
    """The 8/3 lemma covers intervals of length <= 1 only; a bin [1, 2.5] must be rejected."""
    reg, ass, _ = loaded
    cert = json.loads((root / "certificates" / "pc-provencap-A2.json").read_text())
    errs = check_registry._cap_scope_errors("x", cert, ["RH", "lemma-intF-8-3"], ass)
    assert errs == []
    cert["outside"]["edges"] = ["1", "5/2"]
    assert check_registry._cap_scope_errors("x", cert, ["RH", "lemma-intF-8-3"], ass)
    cert["outside"]["edges"], cert["outside"]["caps"] = ["1", "2"], ["2"]   # cap below 8/3
    assert check_registry._cap_scope_errors("x", cert, ["RH", "lemma-intF-8-3"], ass)
    cert["outside"]["caps"] = ["8/3"]
    assert check_registry._cap_scope_errors("x", cert, ["RH"], ass)          # no cited cap assumption


def test_unproven_cap_scale_must_match_the_certificate(loaded, root):
    reg, ass, _ = loaded
    cert = json.loads((root / "certificates" / "pc-legacy-caps-A2-scale4.json").read_text())
    assert check_registry._cap_scope_errors("x", cert, ["caps-scale-4"], ass) == []
    assert check_registry._cap_scope_errors("x", cert, ["caps-scale-8"], ass)   # caps too small for scale 8


def test_every_certificate_records_its_provenance(root):
    for p in (root / "certificates").glob("*.json"):
        prov = json.loads(p.read_text())["provenance"]
        assert check_registry.PROVENANCE_FIELDS <= set(prov)


def test_failed_experiments_and_defects_are_kept(loaded):
    reg, _, _ = loaded
    kinds = [e["kind"] for e in reg["results"]]
    assert kinds.count("negative-result") >= 4 and kinds.count("known-defect") >= 3


def test_research_workflow_does_not_mask_audit_script_failures(root):
    """2026-11-23 independent audit, finding F4: only the documented check3_exact failure may be tolerated."""
    text = (root / ".github" / "workflows" / "research.yml").read_text()
    assert "|| echo" not in text
    assert 'exit "$status"' in text and "SolverError" in text
