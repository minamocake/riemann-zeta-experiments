"""Verify the independent historical v1/v2 archives and their source provenance."""

import hashlib
import json
from pathlib import Path


def test_historical_snapshots_are_complete_and_accounted_for(root: Path):
    manifest = json.loads((root / "legacy/prototypes/SOURCE_MANIFEST.json").read_text())
    for version, expected_count in (("v1_positive_boxes", 13), ("v2_signed_boxes", 27)):
        entry = manifest["archives"][version]
        base = root / "legacy/prototypes" / version
        assert entry["file_count"] == expected_count
        assert len(entry["files"]) == expected_count
        assert all(len(rec["source_sha256"]) == 64 for rec in entry["files"])
        for rec in entry["files"]:
            assert hashlib.sha256((base / rec["path"]).read_bytes()).hexdigest() == rec["archived_sha256"]
            if rec["transformation"] == "byte-identical":
                assert rec["archived_sha256"] == rec["source_sha256"]
        extras = {p.relative_to(base).as_posix() for p in base.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
        expected = {rec["path"] for rec in entry["files"]} | {"ARCHIVE_STATUS.md"}
        assert extras == expected


def test_only_the_signed_v2_historical_note_is_translated(root: Path):
    entry = json.loads((root / "legacy/prototypes/SOURCE_MANIFEST.json").read_text())["archives"]
    changed = [(version, rec["path"]) for version, d in entry.items()
               for rec in d["files"] if rec["source_sha256"] != rec["archived_sha256"]]
    assert changed == [("v2_signed_boxes", "signed_research_notes.md")]


def test_common_legacy_v1_v2_code_is_byte_identical(root: Path):
    base = root / "legacy/prototypes"
    for name in ("run_experiment.py", "src/pair_sensitivity.py", "requirements.txt",
                 "tests/test_pair_sensitivity.py", "research_notes.md"):
        a = (base / "v1_positive_boxes" / name).read_bytes()
        b = (base / "v2_signed_boxes" / name).read_bytes()
        assert a == b, name
