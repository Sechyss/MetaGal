"""
Visualisation module for 16S metabarcode analysis.

Generates publication-quality figures for alpha diversity, beta diversity
(PCoA ordination), and relative abundance bar charts.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for server/HPC use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

logger = logging.getLogger(__name__)

# Colour palette suitable for colour-blind users (Wong 2011)
CB_PALETTE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]


def plot_alpha_diversity(
    alpha_df: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    metric: str = "shannon",
    output_path: str | None = None,
) -> plt.Figure:
    """
    Box-and-strip plot of an alpha diversity metric grouped by a metadata variable.

    Parameters
    ----------
    alpha_df : pd.DataFrame
        Alpha diversity table (samples × metrics).
    metadata : pd.DataFrame
        Sample metadata; index must match ``alpha_df``.
    group_column : str
        Metadata column defining the groups.
    metric : str
        Diversity metric to plot (default: ``'shannon'``).
    output_path : str or None
        If provided, save the figure to this path.

    Returns
    -------
    plt.Figure
    """
    groups = metadata[group_column].unique()
    data = [
        alpha_df.loc[metadata.index[metadata[group_column] == g], metric].dropna()
        for g in groups
    ]

    fig, ax = plt.subplots(figsize=(max(4, len(groups) * 1.2), 5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.5)

    for patch, colour in zip(bp["boxes"], CB_PALETTE):
        patch.set_facecolor(to_rgba(colour, 0.7))

    # Overlay individual data points
    for i, (d, colour) in enumerate(zip(data, CB_PALETTE), start=1):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(d))
        ax.scatter(i + jitter, d, color=colour, s=30, zorder=3, alpha=0.8)

    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=20, ha="right")
    ax.set_ylabel(metric.replace("_", " ").title())
    ax.set_title(f"Alpha diversity – {metric}")
    fig.tight_layout()

    if output_path:
        os.makedirs(Path(output_path).parent, exist_ok=True)
        fig.savefig(output_path, dpi=150)
        logger.info("Alpha diversity plot saved to %s", output_path)

    return fig


def plot_pcoa(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    output_path: str | None = None,
) -> plt.Figure:
    """
    Principal Coordinates Analysis (PCoA) ordination plot.

    Parameters
    ----------
    distance_matrix : pd.DataFrame
        Square dissimilarity matrix (samples × samples).
    metadata : pd.DataFrame
        Sample metadata; index must match the distance matrix.
    group_column : str
        Metadata column defining the groups (used for colouring points).
    output_path : str or None
        If provided, save the figure to this path.

    Returns
    -------
    plt.Figure
    """
    from scipy.spatial.distance import squareform
    from sklearn.manifold import MDS

    samples = distance_matrix.index.tolist()
    dmat = distance_matrix.loc[samples, samples].values

    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    coords = mds.fit_transform(dmat)

    groups = metadata.loc[samples, group_column]
    unique_groups = groups.unique()
    colour_map = {g: CB_PALETTE[i % len(CB_PALETTE)]
                  for i, g in enumerate(unique_groups)}

    fig, ax = plt.subplots(figsize=(6, 5))
    for group in unique_groups:
        mask = groups == group
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=colour_map[group], label=group, s=60, alpha=0.85,
        )

    ax.set_xlabel("PCoA 1")
    ax.set_ylabel("PCoA 2")
    ax.set_title(f"PCoA – Bray-Curtis dissimilarity\nColoured by {group_column}")
    ax.legend(title=group_column, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()

    if output_path:
        os.makedirs(Path(output_path).parent, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info("PCoA plot saved to %s", output_path)

    return fig


def plot_relative_abundance(
    asv_table: pd.DataFrame,
    taxonomy: pd.DataFrame,
    level: str = "Phylum",
    top_n: int = 10,
    output_path: str | None = None,
) -> plt.Figure:
    """
    Stacked bar chart of relative abundance at a given taxonomic level.

    Parameters
    ----------
    asv_table : pd.DataFrame
        ASV table (ASVs as rows, samples as columns).
    taxonomy : pd.DataFrame
        Taxonomy table with ``Feature ID`` and ``Taxon`` columns.
    level : str
        Taxonomic level to collapse to (default: ``'Phylum'``).
    top_n : int
        Keep only the ``top_n`` most abundant taxa; remainder grouped as
        ``'Other'`` (default: 10).
    output_path : str or None
        If provided, save the figure to this path.

    Returns
    -------
    plt.Figure
    """
    level_index = {
        "Domain": 0, "Kingdom": 0,
        "Phylum": 1, "Class": 2, "Order": 3,
        "Family": 4, "Genus": 5, "Species": 6,
    }

    tax_map = taxonomy.set_index("Feature ID")["Taxon"]
    annotated = asv_table.copy()
    annotated["Taxon"] = annotated.index.map(tax_map).fillna("Unclassified")

    idx = level_index.get(level, 1)

    def _extract_level(taxon: str) -> str:
        parts = [p.strip() for p in taxon.split(";")]
        return parts[idx] if idx < len(parts) else "Unclassified"

    annotated["Level"] = annotated["Taxon"].apply(_extract_level)
    grouped = annotated.drop(columns=["Taxon"]).groupby("Level").sum()

    # Relative abundance
    rel = grouped.div(grouped.sum(axis=0), axis=1) * 100

    # Keep only top N taxa
    top_taxa = rel.sum(axis=1).nlargest(top_n).index
    other = rel.drop(index=top_taxa).sum(axis=0)
    rel = rel.loc[top_taxa]
    rel.loc["Other"] = other

    colours = plt.cm.tab20.colors
    fig, ax = plt.subplots(figsize=(max(6, len(rel.columns) * 0.8), 6))
    rel.T.plot(kind="bar", stacked=True, ax=ax, color=colours, width=0.75)

    ax.set_ylabel("Relative abundance (%)")
    ax.set_xlabel("Sample")
    ax.set_title(f"Relative abundance ({level} level)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.legend(title=level, bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=8)
    fig.tight_layout()

    if output_path:
        os.makedirs(Path(output_path).parent, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info("Relative abundance plot saved to %s", output_path)

    return fig
