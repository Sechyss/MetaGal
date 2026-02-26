#!/usr/bin/env python3
"""
MetaGal – 16S metabarcode analysis pipeline.

Usage
-----
    python pipeline.py --config config/config.yaml

The pipeline executes the following steps in order:

1. Quality control  – FastQC + MultiQC on raw reads
2. Primer trimming  – cutadapt
3. Denoising        – DADA2 (Rscript)
4. Taxonomy         – VSEARCH or QIIME 2 sklearn classifier
5. Diversity        – Alpha (Shannon, Simpson) and beta (Bray-Curtis) metrics
6. Visualisation    – Alpha diversity plots, PCoA, relative abundance charts
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import yaml

from metagal.quality_control import quality_control
from metagal.trimming import trim_samples
from metagal.denoising import run_dada2, load_asv_table
from metagal.taxonomy import classify_vsearch, classify_qiime2, load_taxonomy, merge_taxonomy
from metagal.diversity import alpha_diversity, bray_curtis_matrix, save_diversity_results
from metagal.visualisation import (
    plot_alpha_diversity,
    plot_pcoa,
    plot_relative_abundance,
)

import pandas as pd


def setup_logging(log_level: str = "INFO", log_file: str | None = None) -> None:
    """Configure root logger."""
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(Path(log_file).parent, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )


def load_config(config_path: str) -> dict:
    """
    Load and validate a YAML configuration file.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Parsed configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If the configuration file does not exist.
    ValueError
        If required keys are missing.
    """
    if not Path(config_path).exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    required = {"samples", "output_dir"}
    missing = required - set(config.keys())
    if missing:
        raise ValueError(f"Missing required configuration keys: {missing}")

    return config


def load_samples(samples_config: dict | str) -> dict[str, tuple[str, str]]:
    """
    Build a sample dictionary from the config ``samples`` entry.

    Accepts either:
    - A dictionary mapping sample names to ``{r1: ..., r2: ...}``
    - A path to a TSV file with columns: sample, r1, r2

    Parameters
    ----------
    samples_config : dict or str
        Sample configuration from the YAML config.

    Returns
    -------
    dict[str, tuple[str, str]]
        Mapping of sample name → (R1 path, R2 path).
    """
    if isinstance(samples_config, str):
        # TSV manifest
        manifest = pd.read_csv(samples_config, sep="\t")
        return {
            row["sample"]: (row["r1"], row["r2"])
            for _, row in manifest.iterrows()
        }

    return {
        name: (info["r1"], info["r2"])
        for name, info in samples_config.items()
    }


def run_pipeline(config: dict) -> None:
    """
    Execute the full MetaGal pipeline.

    Parameters
    ----------
    config : dict
        Parsed configuration dictionary (see ``config/config.yaml``).
    """
    log = logging.getLogger("metagal.pipeline")
    output_dir = config["output_dir"]
    threads = config.get("threads", 1)

    # ── 1. Load samples ──────────────────────────────────────────────────────
    samples = load_samples(config["samples"])
    log.info("Loaded %d sample(s)", len(samples))

    # ── 2. Quality control ───────────────────────────────────────────────────
    if config.get("steps", {}).get("quality_control", True):
        log.info("=== Step 1: Quality Control ===")
        all_reads = [f for r1, r2 in samples.values() for f in (r1, r2)]
        quality_control(
            input_files=all_reads,
            output_dir=os.path.join(output_dir, "qc"),
            threads=threads,
        )

    # ── 3. Primer trimming ───────────────────────────────────────────────────
    if config.get("steps", {}).get("trimming", True):
        log.info("=== Step 2: Primer Trimming ===")
        trimming_cfg = config.get("trimming", {})
        trimmed = trim_samples(
            samples=samples,
            output_dir=os.path.join(output_dir, "trimmed"),
            primer_set=trimming_cfg.get("primer_set", "515F_806R"),
            forward_primer=trimming_cfg.get("forward_primer"),
            reverse_primer=trimming_cfg.get("reverse_primer"),
            min_length=trimming_cfg.get("min_length", 100),
            error_rate=trimming_cfg.get("error_rate", 0.1),
            discard_untrimmed=trimming_cfg.get("discard_untrimmed", True),
            threads=threads,
        )
    else:
        trimmed = samples

    # ── 4. Denoising (DADA2) ─────────────────────────────────────────────────
    if config.get("steps", {}).get("denoising", True):
        log.info("=== Step 3: Denoising (DADA2) ===")
        dada2_cfg = config.get("dada2", {})
        asv_table_path, rep_seqs_path = run_dada2(
            samples=trimmed,
            output_dir=os.path.join(output_dir, "dada2"),
            trunc_len_f=dada2_cfg.get("trunc_len_f", 230),
            trunc_len_r=dada2_cfg.get("trunc_len_r", 200),
            threads=threads,
        )
    else:
        asv_table_path = config["asv_table"]
        rep_seqs_path = config["rep_seqs"]

    asv_table = load_asv_table(asv_table_path)

    # ── 5. Taxonomy ──────────────────────────────────────────────────────────
    if config.get("steps", {}).get("taxonomy", True):
        log.info("=== Step 4: Taxonomic Classification ===")
        tax_cfg = config.get("taxonomy", {})
        method = tax_cfg.get("method", "vsearch")

        if method == "vsearch":
            taxonomy_path = classify_vsearch(
                rep_seqs=rep_seqs_path,
                database=tax_cfg["database"],
                output_dir=os.path.join(output_dir, "taxonomy"),
                identity=tax_cfg.get("identity", 0.97),
                threads=threads,
            )
        elif method == "qiime2":
            taxonomy_path = classify_qiime2(
                rep_seqs=rep_seqs_path,
                classifier=tax_cfg["classifier"],
                output_dir=os.path.join(output_dir, "taxonomy"),
                threads=threads,
            )
        else:
            raise ValueError(f"Unknown taxonomy method: '{method}'")

        taxonomy = load_taxonomy(taxonomy_path)
        annotated_table = merge_taxonomy(asv_table, taxonomy)
        annotated_table.to_csv(
            os.path.join(output_dir, "taxonomy", "asv_table_annotated.tsv"),
            sep="\t",
        )
    else:
        taxonomy = None

    # ── 6. Diversity analysis ────────────────────────────────────────────────
    if config.get("steps", {}).get("diversity", True):
        log.info("=== Step 5: Diversity Analysis ===")
        alpha_df = alpha_diversity(asv_table)
        beta_df = bray_curtis_matrix(asv_table)
        diversity_dir = os.path.join(output_dir, "diversity")
        save_diversity_results(alpha_df, beta_df, diversity_dir)

    # ── 7. Visualisation ─────────────────────────────────────────────────────
    if config.get("steps", {}).get("visualisation", True):
        log.info("=== Step 6: Visualisation ===")
        viz_dir = os.path.join(output_dir, "figures")
        metadata_path = config.get("metadata")

        if metadata_path and Path(metadata_path).exists():
            metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
            group_col = config.get("group_column", metadata.columns[0])

            plot_alpha_diversity(
                alpha_df=alpha_df,
                metadata=metadata,
                group_column=group_col,
                output_path=os.path.join(viz_dir, "alpha_diversity.png"),
            )
            plot_pcoa(
                distance_matrix=beta_df,
                metadata=metadata,
                group_column=group_col,
                output_path=os.path.join(viz_dir, "pcoa.png"),
            )

        if taxonomy is not None:
            plot_relative_abundance(
                asv_table=asv_table,
                taxonomy=taxonomy,
                level=config.get("taxonomy_level", "Phylum"),
                top_n=config.get("top_n_taxa", 10),
                output_path=os.path.join(viz_dir, "relative_abundance.png"),
            )

    log.info("=== Pipeline complete. Results in: %s ===", output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MetaGal – 16S metabarcode analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        metavar="FILE",
        help="Path to the YAML configuration file.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        metavar="FILE",
        help="Optional path to write log output to a file.",
    )

    args = parser.parse_args()
    setup_logging(args.log_level, args.log_file)

    config = load_config(args.config)
    run_pipeline(config)


if __name__ == "__main__":
    main()
