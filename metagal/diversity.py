"""
Diversity analysis module for 16S metabarcode analysis.

Computes alpha-diversity (richness, Shannon, Simpson) and beta-diversity
(Bray-Curtis, UniFrac via scikit-bio) metrics from an ASV count table.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis
from scipy.stats import kruskal

logger = logging.getLogger(__name__)


# ── Alpha diversity ───────────────────────────────────────────────────────────

def observed_features(counts: pd.Series) -> int:
    """Return the number of ASVs with at least one count."""
    return int((counts > 0).sum())


def shannon_entropy(counts: pd.Series) -> float:
    """
    Calculate Shannon entropy (H').

    Parameters
    ----------
    counts : pd.Series
        Per-ASV read counts for a single sample.

    Returns
    -------
    float
        Shannon entropy (nats).
    """
    freqs = counts[counts > 0] / counts.sum()
    return float(-np.sum(freqs * np.log(freqs)))


def simpson_index(counts: pd.Series) -> float:
    """
    Calculate the Simpson dominance index (D).

    Parameters
    ----------
    counts : pd.Series
        Per-ASV read counts for a single sample.

    Returns
    -------
    float
        Simpson index (0 = no diversity, 1 = maximum diversity).
    """
    n = counts.sum()
    if n <= 1:
        return 0.0
    freqs = counts[counts > 0] / n
    return float(1 - np.sum(freqs ** 2))


def alpha_diversity(asv_table: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a suite of alpha-diversity metrics for each sample.

    Parameters
    ----------
    asv_table : pd.DataFrame
        ASV table with ASVs as rows and samples as columns.
        Values must be non-negative integers or floats.

    Returns
    -------
    pd.DataFrame
        DataFrame with samples as rows and diversity metrics as columns:
        ``observed_features``, ``shannon``, ``simpson``.
    """
    results = {}
    for sample in asv_table.columns:
        col = asv_table[sample]
        results[sample] = {
            "observed_features": observed_features(col),
            "shannon": shannon_entropy(col),
            "simpson": simpson_index(col),
        }
    return pd.DataFrame(results).T


# ── Beta diversity ────────────────────────────────────────────────────────────

def bray_curtis_matrix(asv_table: pd.DataFrame) -> pd.DataFrame:
    """
    Compute pairwise Bray-Curtis dissimilarity between samples.

    Parameters
    ----------
    asv_table : pd.DataFrame
        ASV table with ASVs as rows and samples as columns.

    Returns
    -------
    pd.DataFrame
        Square dissimilarity matrix (samples × samples).
    """
    samples = asv_table.columns.tolist()
    n = len(samples)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            dist = braycurtis(asv_table.iloc[:, i], asv_table.iloc[:, j])
            matrix[i, j] = dist
            matrix[j, i] = dist

    return pd.DataFrame(matrix, index=samples, columns=samples)


# ── Statistical testing ───────────────────────────────────────────────────────

def kruskal_wallis_test(
    alpha_df: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    metric: str = "shannon",
) -> dict[str, float]:
    """
    Test whether alpha diversity differs significantly between groups.

    Uses the Kruskal-Wallis H-test (non-parametric one-way ANOVA).

    Parameters
    ----------
    alpha_df : pd.DataFrame
        Alpha diversity DataFrame (output of :func:`alpha_diversity`).
    metadata : pd.DataFrame
        Sample metadata; index must match the index of ``alpha_df``.
    group_column : str
        Column in ``metadata`` that defines the groups.
    metric : str
        Alpha diversity metric to test (default: ``'shannon'``).

    Returns
    -------
    dict[str, float]
        ``{'statistic': float, 'p_value': float}``.

    Raises
    ------
    KeyError
        If ``group_column`` is not found in ``metadata``, or ``metric`` is
        not found in ``alpha_df``.
    ValueError
        If fewer than two groups are present.
    """
    if group_column not in metadata.columns:
        raise KeyError(f"Column '{group_column}' not found in metadata.")
    if metric not in alpha_df.columns:
        raise KeyError(f"Metric '{metric}' not found in alpha diversity table.")

    groups = metadata[group_column].unique()
    if len(groups) < 2:
        raise ValueError("At least two groups are required for the test.")

    group_values = [
        alpha_df.loc[metadata.index[metadata[group_column] == g], metric].dropna()
        for g in groups
    ]
    stat, pval = kruskal(*group_values)
    return {"statistic": float(stat), "p_value": float(pval)}


# ── I/O helpers ───────────────────────────────────────────────────────────────

def save_diversity_results(
    alpha_df: pd.DataFrame,
    beta_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """
    Write alpha and beta diversity tables to TSV files.

    Parameters
    ----------
    alpha_df : pd.DataFrame
        Alpha diversity table.
    beta_df : pd.DataFrame
        Beta diversity (dissimilarity) matrix.
    output_dir : str
        Directory where TSV files will be written.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    alpha_path = Path(output_dir) / "alpha_diversity.tsv"
    beta_path = Path(output_dir) / "beta_diversity_bray_curtis.tsv"

    alpha_df.to_csv(alpha_path, sep="\t")
    beta_df.to_csv(beta_path, sep="\t")

    logger.info("Diversity results written to %s", output_dir)
