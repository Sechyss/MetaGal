#!/usr/bin/env python3
"""
Epi2me import helper for MetaGal.

Oxford Nanopore's Epi2me software stores basecalled, demultiplexed reads
in a directory tree organised by barcode (e.g. ``fastq_pass/barcode01/``).
Each barcode directory typically contains multiple FASTQ files that need to
be concatenated into a single per-sample file before the MetaGal pipeline
can be run.

A sample-sheet CSV (``sample_sheet.csv``) maps barcode identifiers to
human-readable sample names.  If no sample-sheet is supplied, barcode
directory names are used as sample names.

Usage
-----
::

    python scripts/epi2me_import.py \\
        --epi2me-dir /path/to/epi2me_output \\
        --sample-sheet /path/to/sample_sheet.csv \\
        --output-dir data/raw \\
        --manifest data/manifest_nanopore.tsv

The script will:

1. Discover all barcode directories under ``--epi2me-dir``
   (``fastq_pass/barcodeXX`` or ``barcodeXX`` at the top level).
2. Concatenate every FASTQ / FASTQ.gz file inside each barcode directory
   into a single ``<sample_name>.fastq.gz`` in ``--output-dir``.
3. Write a TSV manifest (``--manifest``) with columns ``sample`` and
   ``reads`` ready for use with ``config_nanopore.yaml``.

Sample-sheet format
-------------------
A CSV file with at least two columns:

+----------+-------------+
| barcode  | sample_name |
+----------+-------------+
| barcode01| patient_A   |
| barcode02| patient_B   |
+----------+-------------+

Column names are configurable via ``--barcode-col`` and ``--name-col``.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# ── Discovery ─────────────────────────────────────────────────────────────────

def find_barcode_dirs(epi2me_dir: str) -> dict[str, Path]:
    """
    Locate per-barcode FASTQ directories inside an Epi2me output folder.

    Searches for sub-directories whose names match the pattern
    ``barcodeNN`` (case-insensitive) under the given root, checking
    both the root itself and a ``fastq_pass`` sub-directory.

    Parameters
    ----------
    epi2me_dir : str
        Root directory of the Epi2me output.

    Returns
    -------
    dict[str, Path]
        Mapping of barcode name (lower-case, e.g. ``"barcode01"``) →
        absolute Path to the barcode directory.

    Raises
    ------
    FileNotFoundError
        If ``epi2me_dir`` does not exist.
    """
    root = Path(epi2me_dir)
    if not root.exists():
        raise FileNotFoundError(f"Epi2me output directory not found: {epi2me_dir}")

    candidates = [root, root / "fastq_pass"]
    barcode_dirs: dict[str, Path] = {}

    for search_dir in candidates:
        if not search_dir.is_dir():
            continue
        for entry in sorted(search_dir.iterdir()):
            if entry.is_dir() and entry.name.lower().startswith("barcode"):
                key = entry.name.lower()
                if key not in barcode_dirs:
                    barcode_dirs[key] = entry

    return barcode_dirs


def collect_fastq_files(barcode_dir: Path) -> list[Path]:
    """
    Return all FASTQ files (``*.fastq``, ``*.fastq.gz``, ``*.fq``,
    ``*.fq.gz``) inside *barcode_dir*, sorted by name.

    Parameters
    ----------
    barcode_dir : Path
        Path to the barcode directory.

    Returns
    -------
    list[Path]
        Sorted list of FASTQ file paths.
    """
    patterns = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    files: list[Path] = []
    for pattern in patterns:
        files.extend(barcode_dir.glob(pattern))
    return sorted(files)


# ── Concatenation ─────────────────────────────────────────────────────────────

def concatenate_fastq(files: list[Path], output_path: Path) -> None:
    """
    Concatenate one or more FASTQ / FASTQ.gz files into a single gzip file.

    Parameters
    ----------
    files : list[Path]
        Input FASTQ files (may be gzip-compressed or plain).
    output_path : Path
        Path of the output ``.fastq.gz`` file.

    Raises
    ------
    ValueError
        If *files* is empty.
    """
    if not files:
        raise ValueError("No FASTQ files to concatenate.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(output_path, "wb") as out_fh:
        for fq in files:
            if fq.name.endswith((".fastq.gz", ".fq.gz")):
                with gzip.open(fq, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)
            else:
                with open(fq, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)

    logger.debug("Concatenated %d file(s) → %s", len(files), output_path)


# ── Sample-sheet parsing ──────────────────────────────────────────────────────

def load_sample_sheet(
    sample_sheet_path: str,
    barcode_col: str = "barcode",
    name_col: str = "sample_name",
) -> dict[str, str]:
    """
    Parse an Epi2me sample-sheet CSV into a barcode → sample-name mapping.

    Parameters
    ----------
    sample_sheet_path : str
        Path to the sample-sheet CSV file.
    barcode_col : str
        Name of the column containing barcode identifiers (default:
        ``"barcode"``).
    name_col : str
        Name of the column containing sample names (default:
        ``"sample_name"``).

    Returns
    -------
    dict[str, str]
        Mapping of lower-case barcode string → sample name.

    Raises
    ------
    FileNotFoundError
        If the sample-sheet file does not exist.
    KeyError
        If the expected columns are missing.
    """
    path = Path(sample_sheet_path)
    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {sample_sheet_path}")

    df = pd.read_csv(path)

    missing = {barcode_col, name_col} - set(df.columns)
    if missing:
        raise KeyError(
            f"Sample sheet is missing column(s): {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    return {
        row[barcode_col].lower(): row[name_col]
        for _, row in df.iterrows()
    }


# ── Main import logic ─────────────────────────────────────────────────────────

def import_epi2me(
    epi2me_dir: str,
    output_dir: str,
    manifest_path: str,
    sample_sheet_path: str | None = None,
    barcode_col: str = "barcode",
    name_col: str = "sample_name",
) -> pd.DataFrame:
    """
    Import Epi2me output into the MetaGal pipeline format.

    Concatenates per-barcode FASTQ files and writes a sample manifest TSV.

    Parameters
    ----------
    epi2me_dir : str
        Root directory of the Epi2me output.
    output_dir : str
        Directory where concatenated FASTQ files will be written.
    manifest_path : str
        Path to write the TSV manifest (``sample``, ``reads`` columns).
    sample_sheet_path : str or None
        Path to the Epi2me sample-sheet CSV.  If ``None``, barcode directory
        names are used as sample names.
    barcode_col : str
        Sample-sheet column for barcode IDs (default: ``"barcode"``).
    name_col : str
        Sample-sheet column for sample names (default: ``"sample_name"``).

    Returns
    -------
    pd.DataFrame
        The written manifest as a DataFrame with columns ``sample`` and
        ``reads``.

    Raises
    ------
    FileNotFoundError
        If ``epi2me_dir`` or ``sample_sheet_path`` does not exist.
    RuntimeError
        If no barcode directories are found.
    """
    logger.info("Scanning Epi2me output directory: %s", epi2me_dir)
    barcode_dirs = find_barcode_dirs(epi2me_dir)
    if not barcode_dirs:
        raise RuntimeError(
            f"No barcode directories found in {epi2me_dir}. "
            "Expected sub-directories named 'barcodeNN' (e.g. barcode01)."
        )
    logger.info("Found %d barcode directory(ies)", len(barcode_dirs))

    # Build barcode → sample-name mapping
    if sample_sheet_path:
        barcode_to_name = load_sample_sheet(sample_sheet_path, barcode_col, name_col)
    else:
        barcode_to_name = {bc: bc for bc in barcode_dirs}

    os.makedirs(output_dir, exist_ok=True)
    records: list[dict[str, str]] = []

    for barcode, bdir in sorted(barcode_dirs.items()):
        sample_name = barcode_to_name.get(barcode, barcode)
        fastq_files = collect_fastq_files(bdir)

        if not fastq_files:
            logger.warning("No FASTQ files found in %s – skipping", bdir)
            continue

        out_fastq = Path(output_dir) / f"{sample_name}.fastq.gz"
        logger.info(
            "Concatenating %d file(s) for '%s' → %s",
            len(fastq_files), sample_name, out_fastq,
        )
        concatenate_fastq(fastq_files, out_fastq)
        records.append({"sample": sample_name, "reads": str(out_fastq)})

    manifest_df = pd.DataFrame(records, columns=["sample", "reads"])
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    logger.info(
        "Manifest written to %s (%d sample(s))", manifest_path, len(manifest_df)
    )
    return manifest_df


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Import Epi2me output into the MetaGal pipeline format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--epi2me-dir", "-e",
        required=True,
        metavar="DIR",
        help="Root directory of the Epi2me output.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        metavar="DIR",
        help="Directory where concatenated FASTQ files will be written.",
    )
    parser.add_argument(
        "--manifest", "-m",
        required=True,
        metavar="FILE",
        help="Path to write the TSV manifest (columns: sample, reads).",
    )
    parser.add_argument(
        "--sample-sheet", "-s",
        default=None,
        metavar="FILE",
        help=(
            "Epi2me sample-sheet CSV mapping barcodes to sample names. "
            "If omitted, barcode directory names are used."
        ),
    )
    parser.add_argument(
        "--barcode-col",
        default="barcode",
        metavar="COL",
        help="Sample-sheet column for barcode IDs (default: barcode).",
    )
    parser.add_argument(
        "--name-col",
        default="sample_name",
        metavar="COL",
        help="Sample-sheet column for sample names (default: sample_name).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    import_epi2me(
        epi2me_dir=args.epi2me_dir,
        output_dir=args.output_dir,
        manifest_path=args.manifest,
        sample_sheet_path=args.sample_sheet,
        barcode_col=args.barcode_col,
        name_col=args.name_col,
    )


if __name__ == "__main__":
    main()
