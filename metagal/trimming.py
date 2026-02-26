"""
Primer trimming module for 16S metabarcode analysis.

Uses cutadapt to remove primer sequences from paired-end reads.
Supports standard 16S primer sets (515F/806R, 27F/338R, etc.).
"""

import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Common 16S primer sequences
PRIMERS = {
    "515F_806R": {
        "forward": "GTGYCAGCMGCCGCGGTAA",
        "reverse": "GGACTACNVGGGTWTCTAAT",
        "region": "V4",
    },
    "27F_338R": {
        "forward": "AGAGTTTGATCMTGGCTCAG",
        "reverse": "TGCTGCCTCCCGTAGGAGT",
        "region": "V1-V2",
    },
    "341F_806R": {
        "forward": "CCTACGGGNGGCWGCAG",
        "reverse": "GGACTACNVGGGTWTCTAAT",
        "region": "V3-V4",
    },
}


def trim_primers(
    forward_reads: str,
    reverse_reads: str,
    output_dir: str,
    forward_primer: str,
    reverse_primer: str,
    min_length: int = 100,
    error_rate: float = 0.1,
    discard_untrimmed: bool = True,
    threads: int = 1,
) -> tuple[str, str]:
    """
    Trim primer sequences from paired-end reads using cutadapt.

    Parameters
    ----------
    forward_reads : str
        Path to forward (R1) FASTQ file.
    reverse_reads : str
        Path to reverse (R2) FASTQ file.
    output_dir : str
        Directory where trimmed reads will be written.
    forward_primer : str
        Forward primer nucleotide sequence (5'→3').
    reverse_primer : str
        Reverse primer nucleotide sequence (5'→3').
    min_length : int
        Discard reads shorter than this length after trimming (default: 100).
    error_rate : float
        Maximum allowed error rate for primer matching (default: 0.1).
    discard_untrimmed : bool
        If True, discard read pairs where primers are not found (default: True).
    threads : int
        Number of threads to use (default: 1).

    Returns
    -------
    tuple[str, str]
        Paths to trimmed forward and reverse FASTQ files.

    Raises
    ------
    FileNotFoundError
        If an input file does not exist.
    subprocess.CalledProcessError
        If cutadapt exits with a non-zero return code.
    """
    for f in (forward_reads, reverse_reads):
        if not Path(f).exists():
            raise FileNotFoundError(f"Input file not found: {f}")

    os.makedirs(output_dir, exist_ok=True)

    sample_name = Path(forward_reads).name.replace("_R1", "").replace("_1", "")
    sample_name = sample_name.split(".")[0]

    out_r1 = os.path.join(output_dir, f"{sample_name}_R1_trimmed.fastq.gz")
    out_r2 = os.path.join(output_dir, f"{sample_name}_R2_trimmed.fastq.gz")
    log_file = os.path.join(output_dir, f"{sample_name}_cutadapt.log")

    # Build reverse complement of primers for the paired-end adapter trimming
    rc_forward = _reverse_complement(forward_primer)
    rc_reverse = _reverse_complement(reverse_primer)

    cmd = [
        "cutadapt",
        "-g", forward_primer,
        "-G", reverse_primer,
        "-a", rc_reverse,
        "-A", rc_forward,
        "--minimum-length", str(min_length),
        "--error-rate", str(error_rate),
        "--cores", str(threads),
        "-o", out_r1,
        "-p", out_r2,
    ]

    if discard_untrimmed:
        cmd.append("--discard-untrimmed")

    cmd += [forward_reads, reverse_reads]

    logger.info("Trimming primers from %s", sample_name)
    logger.debug("Command: %s", " ".join(cmd))

    with open(log_file, "w") as log_fh:
        subprocess.run(cmd, check=True, stdout=log_fh, stderr=subprocess.STDOUT)

    logger.info("Trimming complete. Output: %s, %s", out_r1, out_r2)
    return out_r1, out_r2


def trim_samples(
    samples: dict[str, tuple[str, str]],
    output_dir: str,
    primer_set: str = "515F_806R",
    forward_primer: str | None = None,
    reverse_primer: str | None = None,
    min_length: int = 100,
    error_rate: float = 0.1,
    discard_untrimmed: bool = True,
    threads: int = 1,
) -> dict[str, tuple[str, str]]:
    """
    Trim primers from multiple samples.

    Parameters
    ----------
    samples : dict[str, tuple[str, str]]
        Dictionary mapping sample names to (R1, R2) file path tuples.
    output_dir : str
        Directory where trimmed reads will be written.
    primer_set : str
        Name of a built-in primer set to use (default: '515F_806R').
        Ignored if ``forward_primer`` and ``reverse_primer`` are provided.
    forward_primer : str or None
        Custom forward primer sequence. Overrides ``primer_set`` if provided.
    reverse_primer : str or None
        Custom reverse primer sequence. Overrides ``primer_set`` if provided.
    min_length : int
        Minimum read length after trimming (default: 100).
    error_rate : float
        Maximum error rate for primer matching (default: 0.1).
    discard_untrimmed : bool
        Discard read pairs without detected primers (default: True).
    threads : int
        Number of threads per cutadapt call (default: 1).

    Returns
    -------
    dict[str, tuple[str, str]]
        Dictionary mapping sample names to trimmed (R1, R2) file paths.

    Raises
    ------
    ValueError
        If ``primer_set`` is not recognised and custom primers are not provided.
    """
    if forward_primer is None or reverse_primer is None:
        if primer_set not in PRIMERS:
            raise ValueError(
                f"Unknown primer set '{primer_set}'. "
                f"Choose from: {list(PRIMERS.keys())} "
                "or supply forward_primer and reverse_primer."
            )
        fwd = PRIMERS[primer_set]["forward"]
        rev = PRIMERS[primer_set]["reverse"]
    else:
        fwd = forward_primer
        rev = reverse_primer

    trimmed = {}
    for sample, (r1, r2) in samples.items():
        sample_out = os.path.join(output_dir, sample)
        trimmed[sample] = trim_primers(
            r1, r2, sample_out, fwd, rev,
            min_length=min_length,
            error_rate=error_rate,
            discard_untrimmed=discard_untrimmed,
            threads=threads,
        )

    return trimmed


# ── Helpers ──────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgtNnYyRrSsWwKkMmBbDdHhVv",
                             "TGCAtgcaNnRrYySsWwMmKkVvHhDdBb")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (IUPAC bases supported)."""
    return seq.translate(_COMPLEMENT)[::-1]
