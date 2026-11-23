"""Fail if any tracked text file or path contains Japanese, Chinese or Korean script.

Mathematical symbols and accented Latin letters (e.g. in author names) are allowed.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CJK = re.compile("[\\u3040-\\u30ff\\u3400-\\u4dbf\\u4e00-\\u9fff\\uac00-\\ud7af\\uff66-\\uff9f\\u3000-\\u303f]")
# Generated, untracked directories. .ruff_cache holds binary cache files without a suffix; decoding
# them as UTF-8 can produce CJK code points and false failures (CI on main, 1e4d6d1).
SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".ruff_cache", "outputs", ".venv"}
TEXT_SUFFIXES = {".py", ".md", ".yaml", ".yml", ".json", ".toml", ".cff", ".txt", ".csv", ".cfg", ".sha256", ""}


def find_problems(root: Path = ROOT) -> list[str]:
    problems = []
    for p in root.rglob("*"):
        rel = p.relative_to(root)
        if set(rel.parts) & SKIP_DIRS:
            continue
        if CJK.search(str(rel)):
            problems.append(f"non-English path: {rel}")
        if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES:
            for n, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if CJK.search(line):
                    problems.append(f"{rel}:{n}: non-English text")
    return problems


def main(root: Path = ROOT) -> int:
    problems = find_problems(root)
    for msg in problems:
        print(msg)
    print("English-only check " + ("FAILED" if problems else "passed"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
