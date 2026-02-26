"""
Quality control module for 16S metabarcode analysis.

Runs FastQC and MultiQC on raw FASTQ files to assess sequencing quality
before downstream processing.
"""

import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_fastqc(input_files: list[str], output_dir: str, threads: int = 1) -> None:
    """
    Run FastQC on a list of FASTQ files.

    Parameters
    ----------
    input_files : list[str]
        Paths to input FASTQ files.
    output_dir : str
        Directory where FastQC reports will be written.
    threads : int
        Number of threads to use (default: 1).

    Raises
    ------
    FileNotFoundError
        If an input file does not exist.
    subprocess.CalledProcessError
        If FastQC exits with a non-zero return code.
    """
    os.makedirs(output_dir, exist_ok=True)

    for f in input_files:
        if not Path(f).exists():
            raise FileNotFoundError(f"Input file not found: {f}")

    cmd = [
        "fastqc",
        "--outdir", output_dir,
        "--threads", str(threads),
    ] + input_files

    logger.info("Running FastQC on %d file(s)", len(input_files))
    logger.debug("Command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("FastQC complete. Reports written to %s", output_dir)


def run_multiqc(input_dir: str, output_dir: str) -> None:
    """
    Run MultiQC to aggregate FastQC reports.

    Parameters
    ----------
    input_dir : str
        Directory containing FastQC reports.
    output_dir : str
        Directory where the MultiQC report will be written.

    Raises
    ------
    subprocess.CalledProcessError
        If MultiQC exits with a non-zero return code.
    """
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "multiqc",
        input_dir,
        "--outdir", output_dir,
        "--force",
    ]

    logger.info("Running MultiQC")
    logger.debug("Command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("MultiQC complete. Report written to %s", output_dir)


def quality_control(
    input_files: list[str],
    output_dir: str,
    threads: int = 1,
) -> str:
    """
    Run full quality control workflow (FastQC + MultiQC).

    Parameters
    ----------
    input_files : list[str]
        Paths to raw FASTQ files.
    output_dir : str
        Root output directory for QC reports.
    threads : int
        Number of threads to use for FastQC (default: 1).

    Returns
    -------
    str
        Path to the MultiQC report directory.
    """
    fastqc_dir = os.path.join(output_dir, "fastqc")
    multiqc_dir = os.path.join(output_dir, "multiqc")

    run_fastqc(input_files, fastqc_dir, threads=threads)
    run_multiqc(fastqc_dir, multiqc_dir)

    return multiqc_dir
