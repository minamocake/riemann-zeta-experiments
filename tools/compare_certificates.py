"""Compare the mathematical content of two certificate directories (provenance is ignored).

Usage: python tools/compare_certificates.py COMMITTED_DIR REGENERATED_DIR
Exit status 1 if any certificate is missing or differs in widths, components, outside or claim.
"""
import json
import sys
from pathlib import Path

KEYS = ("widths", "components", "outside", "claimed")


def main(a: str, b: str) -> int:
    bad = 0
    for p in sorted(Path(a).glob("*.json")):
        q = Path(b) / p.name
        if not q.is_file():
            print(f"missing after regeneration: {p.name}")
            bad += 1
            continue
        x, y = json.loads(p.read_text()), json.loads(q.read_text())
        diff = [k for k in KEYS if x.get(k) != y.get(k)]
        print(f"{p.name:36s} " + ("reproduced" if not diff else f"DIFFERS in {diff}"))
        bad += bool(diff)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
