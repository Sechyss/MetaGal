"""
Taxonomic classification module for 16S metabarcode analysis.

Assigns taxonomy to representative ASV sequences using a reference database
(SILVA, GreenGenes 2, or RDP) via VSEARCH or the QIIME 2 sklearn classifier.
"""

import logging
import os
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Recognised classification methods
METHODS = ("vsearch", "qiime2")


def classify_vsearch(
    rep_seqs: str,
    database: str,
    output_dir: str,
    identity: float = 0.97,
    threads: int = 1,
) -> str:
    """
    Assign taxonomy using VSEARCH global-alignment search.

    Parameters
    ----------
    rep_seqs : str
        Path to representative sequences FASTA (e.g. output of DADA2).
    database : str
        Path to the reference FASTA database with taxonomy in headers
        (SILVA or GreenGenes format).
    output_dir : str
        Directory where taxonomy results will be written.
    identity : float
        Minimum identity threshold for a hit to be accepted (default: 0.97).
    threads : int
        Number of threads (default: 1).

    Returns
    -------
    str
        Path to the taxonomy TSV file.

    Raises
    ------
    FileNotFoundError
        If the representative sequences or database file do not exist.
    subprocess.CalledProcessError
        If VSEARCH exits with a non-zero return code.
    """
    for f in (rep_seqs, database):
        if not Path(f).exists():
            raise FileNotFoundError(f"File not found: {f}")

    os.makedirs(output_dir, exist_ok=True)
    taxonomy_path = os.path.join(output_dir, "taxonomy.tsv")
    no_hit_path = os.path.join(output_dir, "no_hits.fasta")

    cmd = [
        "vsearch",
        "--usearch_global", rep_seqs,
        "--db", database,
        "--id", str(identity),
        "--blast6out", taxonomy_path,
        "--notmatched", no_hit_path,
        "--threads", str(threads),
        "--top_hits_only",
    ]

    logger.info("Running VSEARCH taxonomy classification")
    logger.debug("Command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("Taxonomy written to %s", taxonomy_path)
    return taxonomy_path


def classify_qiime2(
    rep_seqs: str,
    classifier: str,
    output_dir: str,
    threads: int = 1,
) -> str:
    """
    Assign taxonomy using a pre-trained QIIME 2 Naive-Bayes classifier.

    Parameters
    ----------
    rep_seqs : str
        Path to representative sequences FASTA.
    classifier : str
        Path to the QIIME 2 classifier artefact (.qza).
    output_dir : str
        Directory where taxonomy results will be written.
    threads : int
        Number of threads (default: 1).

    Returns
    -------
    str
        Path to the taxonomy TSV file exported from QIIME 2.

    Raises
    ------
    FileNotFoundError
        If the representative sequences or classifier file do not exist.
    subprocess.CalledProcessError
        If a QIIME 2 command exits with a non-zero return code.
    """
    for f in (rep_seqs, classifier):
        if not Path(f).exists():
            raise FileNotFoundError(f"File not found: {f}")

    os.makedirs(output_dir, exist_ok=True)

    # Import representative sequences into QIIME 2 artefact
    rep_seqs_qza = os.path.join(output_dir, "rep_seqs.qza")
    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", rep_seqs,
        "--output-path", rep_seqs_qza,
    ], check=True)

    # Run classifier
    taxonomy_qza = os.path.join(output_dir, "taxonomy.qza")
    subprocess.run([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", classifier,
        "--i-reads", rep_seqs_qza,
        "--o-classification", taxonomy_qza,
        "--p-n-jobs", str(threads),
    ], check=True)

    # Export taxonomy to TSV
    export_dir = os.path.join(output_dir, "taxonomy_export")
    subprocess.run([
        "qiime", "tools", "export",
        "--input-path", taxonomy_qza,
        "--output-path", export_dir,
    ], check=True)

    taxonomy_path = os.path.join(export_dir, "taxonomy.tsv")
    logger.info("Taxonomy written to %s", taxonomy_path)
    return taxonomy_path


def load_taxonomy(taxonomy_path: str) -> pd.DataFrame:
    """
    Load a taxonomy file into a DataFrame.

    Supports VSEARCH BLAST-6 format and QIIME 2 TSV format.

    Parameters
    ----------
    taxonomy_path : str
        Path to the taxonomy TSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with at least ``Feature ID`` and ``Taxon`` columns.

    Raises
    ------
    FileNotFoundError
        If the taxonomy file does not exist.
    """
    if not Path(taxonomy_path).exists():
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_path}")

    # Detect whether the file has a QIIME 2 header row (Feature ID\tTaxon\t...)
    # or is a headerless VSEARCH BLAST-6 file.
    with open(taxonomy_path) as fh:
        first_line = fh.readline()
    has_header = first_line.startswith("Feature ID")

    raw = pd.read_csv(
        taxonomy_path,
        sep="\t",
        header=0 if has_header else None,
        usecols=[0, 1],
    )
    raw.columns = ["Feature ID", "Taxon"]
    return raw


def merge_taxonomy(
    asv_table: pd.DataFrame,
    taxonomy: pd.DataFrame,
) -> pd.DataFrame:
    """
    Annotate an ASV table with taxonomic assignments.

    Parameters
    ----------
    asv_table : pd.DataFrame
        ASV table with ASV identifiers as the index.
    taxonomy : pd.DataFrame
        Taxonomy DataFrame with ``Feature ID`` and ``Taxon`` columns.

    Returns
    -------
    pd.DataFrame
        ASV table with an additional ``Taxon`` column.
    """
    tax_map = taxonomy.set_index("Feature ID")["Taxon"]
    annotated = asv_table.copy()
    annotated.insert(0, "Taxon", annotated.index.map(tax_map))
    return annotated
