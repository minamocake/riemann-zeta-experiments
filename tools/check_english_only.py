"""Fail if any tracked text file or path contains Japanese, Chinese or Korean script.

Mathematical symbols and accented Latin letters (e.g. in author names) are allowed.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CJK = re.compile("[\\u3040-\\u30ff\\u3400-\\u4dbf\\u4e00-\\u9fff\\uac00-\\ud7af\\uff66-\\uff9f\\u3000-\\u303f]")
SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", "outputs", ".venv"}
TEXT_SUFFIXES = {".py", ".md", ".yaml", ".yml", ".json", ".toml", ".cff", ".txt", ".csv", ".cfg", ".sha256", ""}


def main() -> int:
    problems = []
    for p in ROOT.rglob("*"):
        rel = p.relative_to(ROOT)
        if set(rel.parts) & SKIP_DIRS:
            continue
        if CJK.search(str(rel)):
            problems.append(f"non-English path: {rel}")
        if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES:
            for n, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if CJK.search(line):
                    problems.append(f"{rel}:{n}: non-English text")
    for msg in problems:
        print(msg)
    print("English-only check " + ("FAILED" if problems else "passed"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
