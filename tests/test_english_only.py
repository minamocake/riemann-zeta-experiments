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
