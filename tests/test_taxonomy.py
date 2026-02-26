"""
Tests for metagal.taxonomy – taxonomy table loading and ASV annotation.
"""

import io

import pandas as pd
import pytest

from metagal.taxonomy import load_taxonomy, merge_taxonomy


# ── load_taxonomy ─────────────────────────────────────────────────────────────

def test_load_taxonomy_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_taxonomy(str(tmp_path / "nonexistent.tsv"))


def test_load_taxonomy_blast6(tmp_path):
    """VSEARCH BLAST-6 output: two tab-separated columns (Feature ID, Taxon)."""
    tsv = tmp_path / "taxonomy.tsv"
    tsv.write_text(
        "ASV1\tBacteria;Proteobacteria;Gammaproteobacteria\n"
        "ASV2\tBacteria;Firmicutes;Bacilli\n"
    )
    df = load_taxonomy(str(tsv))
    assert "Feature ID" in df.columns
    assert "Taxon" in df.columns
    assert len(df) == 2
    assert df.iloc[0]["Feature ID"] == "ASV1"


# ── merge_taxonomy ────────────────────────────────────────────────────────────

@pytest.fixture()
def asv_table():
    return pd.DataFrame(
        {"S1": [10, 5], "S2": [3, 7]},
        index=["ASV1", "ASV2"],
    )


@pytest.fixture()
def taxonomy_df():
    return pd.DataFrame({
        "Feature ID": ["ASV1", "ASV2"],
        "Taxon": [
            "Bacteria;Proteobacteria;Gammaproteobacteria",
            "Bacteria;Firmicutes;Bacilli",
        ],
    })


def test_merge_taxonomy_adds_column(asv_table, taxonomy_df):
    result = merge_taxonomy(asv_table, taxonomy_df)
    assert "Taxon" in result.columns


def test_merge_taxonomy_values(asv_table, taxonomy_df):
    result = merge_taxonomy(asv_table, taxonomy_df)
    assert result.loc["ASV1", "Taxon"] == "Bacteria;Proteobacteria;Gammaproteobacteria"
    assert result.loc["ASV2", "Taxon"] == "Bacteria;Firmicutes;Bacilli"


def test_merge_taxonomy_preserves_counts(asv_table, taxonomy_df):
    result = merge_taxonomy(asv_table, taxonomy_df)
    assert result.loc["ASV1", "S1"] == 10
    assert result.loc["ASV2", "S2"] == 7


def test_merge_taxonomy_unknown_asv(asv_table):
    """ASVs without a taxonomy entry should have NaN."""
    tax = pd.DataFrame({
        "Feature ID": ["ASV1"],
        "Taxon": ["Bacteria;Proteobacteria"],
    })
    result = merge_taxonomy(asv_table, tax)
    assert pd.isna(result.loc["ASV2", "Taxon"])
