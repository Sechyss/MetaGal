"""
Tests for metagal.trimming – primer handling utilities.
"""

import pytest

from metagal.trimming import _reverse_complement, PRIMERS


# ── _reverse_complement ───────────────────────────────────────────────────────

@pytest.mark.parametrize("seq,expected", [
    ("ATCG", "CGAT"),
    ("AAAA", "TTTT"),
    ("GCGC", "GCGC"),
    ("A",    "T"),
    ("",     ""),
    ("ATCGN", "NCGAT"),
])
def test_reverse_complement(seq, expected):
    assert _reverse_complement(seq) == expected


def test_reverse_complement_iupac():
    """Degenerate IUPAC bases should be complemented correctly."""
    # Y (C/T) → R (A/G), R (A/G) → Y (C/T)
    assert _reverse_complement("YR") == "YR"
    # S (G/C) → S (G/C)
    assert _reverse_complement("S") == "S"


# ── PRIMERS dictionary ────────────────────────────────────────────────────────

def test_primers_keys():
    for key, info in PRIMERS.items():
        assert "forward" in info
        assert "reverse" in info
        assert "region" in info


def test_primers_non_empty_sequences():
    for key, info in PRIMERS.items():
        assert len(info["forward"]) > 0
        assert len(info["reverse"]) > 0
