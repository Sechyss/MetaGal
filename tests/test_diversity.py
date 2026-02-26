"""
Tests for metagal.diversity – alpha and beta diversity calculations.
"""

import numpy as np
import pandas as pd
import pytest

from metagal.diversity import (
    observed_features,
    shannon_entropy,
    simpson_index,
    alpha_diversity,
    bray_curtis_matrix,
    kruskal_wallis_test,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture()
def simple_asv_table():
    """3 ASVs × 4 samples."""
    return pd.DataFrame(
        {
            "S1": [10, 0, 0],
            "S2": [5, 5, 0],
            "S3": [3, 3, 4],
            "S4": [0, 0, 10],
        },
        index=["ASV1", "ASV2", "ASV3"],
    )


@pytest.fixture()
def metadata():
    return pd.DataFrame(
        {"treatment": ["A", "A", "B", "B"]},
        index=["S1", "S2", "S3", "S4"],
    )


# ── observed_features ─────────────────────────────────────────────────────────

def test_observed_features_all_present():
    counts = pd.Series([1, 2, 3])
    assert observed_features(counts) == 3


def test_observed_features_some_absent():
    counts = pd.Series([0, 5, 0, 3])
    assert observed_features(counts) == 2


def test_observed_features_all_absent():
    counts = pd.Series([0, 0, 0])
    assert observed_features(counts) == 0


# ── shannon_entropy ───────────────────────────────────────────────────────────

def test_shannon_single_species():
    """One taxon → zero diversity."""
    counts = pd.Series([100, 0, 0])
    assert shannon_entropy(counts) == pytest.approx(0.0)


def test_shannon_equal_abundance():
    """Two equally abundant taxa → ln(2)."""
    counts = pd.Series([50, 50])
    assert shannon_entropy(counts) == pytest.approx(np.log(2), rel=1e-6)


def test_shannon_three_equal():
    counts = pd.Series([10, 10, 10])
    assert shannon_entropy(counts) == pytest.approx(np.log(3), rel=1e-6)


# ── simpson_index ─────────────────────────────────────────────────────────────

def test_simpson_single_species():
    counts = pd.Series([100, 0])
    assert simpson_index(counts) == pytest.approx(0.0, abs=1e-9)


def test_simpson_equal_two():
    counts = pd.Series([1, 1])
    assert simpson_index(counts) == pytest.approx(0.5, rel=1e-6)


def test_simpson_range():
    for n in range(2, 10):
        val = simpson_index(pd.Series([1] * n))
        assert 0.0 <= val <= 1.0


def test_simpson_empty():
    assert simpson_index(pd.Series([0])) == 0.0


# ── alpha_diversity ───────────────────────────────────────────────────────────

def test_alpha_diversity_shape(simple_asv_table):
    result = alpha_diversity(simple_asv_table)
    assert result.shape == (4, 3)
    assert set(result.columns) == {"observed_features", "shannon", "simpson"}


def test_alpha_diversity_monoculture(simple_asv_table):
    """S1 has only one ASV → zero Shannon and Simpson."""
    result = alpha_diversity(simple_asv_table)
    assert result.loc["S1", "shannon"] == pytest.approx(0.0, abs=1e-9)
    assert result.loc["S1", "simpson"] == pytest.approx(0.0, abs=1e-9)


def test_alpha_diversity_max_richness(simple_asv_table):
    """S3 has the most even distribution → highest Shannon."""
    result = alpha_diversity(simple_asv_table)
    assert result.loc["S3", "shannon"] >= result.loc["S1", "shannon"]
    assert result.loc["S3", "shannon"] >= result.loc["S2", "shannon"]


# ── bray_curtis_matrix ────────────────────────────────────────────────────────

def test_bray_curtis_shape(simple_asv_table):
    matrix = bray_curtis_matrix(simple_asv_table)
    assert matrix.shape == (4, 4)


def test_bray_curtis_diagonal_zero(simple_asv_table):
    matrix = bray_curtis_matrix(simple_asv_table)
    for sample in matrix.index:
        assert matrix.loc[sample, sample] == pytest.approx(0.0)


def test_bray_curtis_symmetry(simple_asv_table):
    matrix = bray_curtis_matrix(simple_asv_table)
    for i in matrix.index:
        for j in matrix.columns:
            assert matrix.loc[i, j] == pytest.approx(matrix.loc[j, i])


def test_bray_curtis_range(simple_asv_table):
    matrix = bray_curtis_matrix(simple_asv_table)
    assert (matrix.values >= 0).all()
    assert (matrix.values <= 1).all()


def test_bray_curtis_identical_samples():
    """Two identical samples have distance 0."""
    table = pd.DataFrame({"A": [1, 2, 3], "B": [1, 2, 3]}, index=["X", "Y", "Z"])
    matrix = bray_curtis_matrix(table)
    assert matrix.loc["A", "B"] == pytest.approx(0.0)


def test_bray_curtis_disjoint_samples():
    """Two samples with no shared taxa have distance 1."""
    table = pd.DataFrame({"A": [5, 0], "B": [0, 5]}, index=["X", "Y"])
    matrix = bray_curtis_matrix(table)
    assert matrix.loc["A", "B"] == pytest.approx(1.0)


# ── kruskal_wallis_test ───────────────────────────────────────────────────────

def test_kruskal_wallis_returns_keys(simple_asv_table, metadata):
    alpha_df = alpha_diversity(simple_asv_table)
    result = kruskal_wallis_test(alpha_df, metadata, "treatment", "shannon")
    assert "statistic" in result
    assert "p_value" in result


def test_kruskal_wallis_missing_group_column(simple_asv_table, metadata):
    alpha_df = alpha_diversity(simple_asv_table)
    with pytest.raises(KeyError):
        kruskal_wallis_test(alpha_df, metadata, "nonexistent", "shannon")


def test_kruskal_wallis_missing_metric(simple_asv_table, metadata):
    alpha_df = alpha_diversity(simple_asv_table)
    with pytest.raises(KeyError):
        kruskal_wallis_test(alpha_df, metadata, "treatment", "nonexistent")


def test_kruskal_wallis_single_group(simple_asv_table):
    meta = pd.DataFrame(
        {"treatment": ["A", "A", "A", "A"]},
        index=["S1", "S2", "S3", "S4"],
    )
    alpha_df = alpha_diversity(simple_asv_table)
    with pytest.raises(ValueError):
        kruskal_wallis_test(alpha_df, meta, "treatment", "shannon")
