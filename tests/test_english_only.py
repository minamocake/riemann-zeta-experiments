"""The repository must be English-only (docs/decisions/0004-english-only.md)."""
import check_english_only


def text(*codepoints):
    return "".join(chr(c) for c in codepoints)


def test_detector_flags_japanese_but_not_mathematical_symbols():
    assert check_english_only.CJK.search(text(0x65E5, 0x672C, 0x8A9E))          # kanji
    assert check_english_only.CJK.search(text(0x3072, 0x3089, 0x304C, 0x306A))  # hiragana
    assert check_english_only.CJK.search(text(0x30AB, 0x30BF, 0x30AB, 0x30CA))  # katakana
    assert not check_english_only.CJK.search(text(0x2264, 0x20, 0x222B, 0x2229, 0x3B1, 0xE7))  # math, Greek, accents


def test_repository_contains_no_non_english_text():
    assert check_english_only.main() == 0


def test_generated_ruff_cache_is_ignored_but_real_text_is_detected(tmp_path):
    """Regression (CI on main, 1e4d6d1): binary Ruff cache files decoded as UTF-8 can contain CJK code
    points. They must be skipped, while non-English text in source files is still reported."""
    cache = tmp_path / ".ruff_cache" / "0.17.6"
    cache.mkdir(parents=True)
    (cache / "16622227928884725061").write_bytes(b"\x00\x01" + text(0x65E5, 0x672C).encode("utf-8") + b"\xff\xfe")
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "clean.py").write_text("x = 1  # English only\n", encoding="utf-8")
    assert check_english_only.find_problems(tmp_path) == []
    (tmp_path / "src" / "bad.py").write_text("# " + text(0x65E5, 0x672C, 0x8A9E) + "\n", encoding="utf-8")
    problems = check_english_only.find_problems(tmp_path)
    assert problems == ["src/bad.py:1: non-English text"]
    assert check_english_only.main(tmp_path) == 1
