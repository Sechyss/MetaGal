"""
Primer trimming module for 16S metabarcode analysis.

Uses cutadapt to remove primer sequences from paired-end Illumina reads.
For Nanopore reads, uses NanoFilt for quality and length filtering, with
an optional cutadapt pass for primer removal.
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


def trim_nanopore_reads(
    reads: str,
    output_dir: str,
    min_length: int = 200,
    max_length: int | None = None,
    min_quality: float = 8.0,
    forward_primer: str | None = None,
    reverse_primer: str | None = None,
    primer_set: str | None = None,
    threads: int = 1,
) -> str:
    """
    Quality-filter and optionally trim primers from Nanopore reads.

    Uses NanoFilt for quality/length filtering and, if primers are
    specified, a subsequent cutadapt pass (single-end) to remove them.

    Parameters
    ----------
    reads : str
        Path to the Nanopore FASTQ file (.fastq or .fastq.gz).
    output_dir : str
        Directory where the filtered reads will be written.
    min_length : int
        Minimum read length to keep (default: 200).
    max_length : int or None
        Maximum read length to keep; ``None`` means no upper limit
        (default: None).
    min_quality : float
        Minimum mean read quality score (Phred) to keep (default: 8.0).
    forward_primer : str or None
        Forward primer sequence for optional primer trimming.
        If *both* ``forward_primer`` and ``reverse_primer`` are ``None``
        and ``primer_set`` is also ``None``, the cutadapt step is skipped.
    reverse_primer : str or None
        Reverse primer sequence for optional primer trimming.
    primer_set : str or None
        Built-in primer set name (see :data:`PRIMERS`).  Overridden by
        explicit ``forward_primer`` / ``reverse_primer`` values.
    threads : int
        Number of threads for cutadapt (default: 1).

    Returns
    -------
    str
        Path to the filtered (and optionally primer-trimmed) FASTQ file.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If ``primer_set`` is specified but not recognised.
    subprocess.CalledProcessError
        If NanoFilt or cutadapt exits with a non-zero return code.
    """
    if not Path(reads).exists():
        raise FileNotFoundError(f"Input file not found: {reads}")

    # Validate primer set early (before running any external tools)
    fwd: str | None = forward_primer
    rev: str | None = reverse_primer

    if fwd is None and rev is None and primer_set is not None:
        if primer_set not in PRIMERS:
            raise ValueError(
                f"Unknown primer set '{primer_set}'. "
                f"Choose from: {list(PRIMERS.keys())} "
                "or supply forward_primer and reverse_primer."
            )
        fwd = PRIMERS[primer_set]["forward"]
        rev = PRIMERS[primer_set]["reverse"]

    os.makedirs(output_dir, exist_ok=True)

    sample_name = Path(reads).name.split(".")[0]
    nanofilt_out = os.path.join(output_dir, f"{sample_name}_nanofilt.fastq.gz")

    # ── NanoFilt: quality + length filtering ──────────────────────────────────
    nanofilt_cmd = [
        "NanoFilt",
        "-q", str(min_quality),
        "-l", str(min_length),
    ]
    if max_length is not None:
        nanofilt_cmd += ["--maxlength", str(max_length)]
    nanofilt_cmd.append(reads)

    logger.info("Running NanoFilt on %s", sample_name)
    logger.debug("NanoFilt command: %s", " ".join(nanofilt_cmd))

    with open(nanofilt_out, "wb") as out_fh:
        nanofilt_proc = subprocess.run(
            nanofilt_cmd,
            check=True,
            stdout=subprocess.PIPE,
        )
        # Compress on the fly via gzip
        gzip_proc = subprocess.run(
            ["gzip", "-c"],
            input=nanofilt_proc.stdout,
            stdout=out_fh,
            check=True,
        )

    output_path = nanofilt_out

    # ── Optional cutadapt primer trimming (single-end) ────────────────────────
    if fwd is not None or rev is not None:
        trimmed_out = os.path.join(output_dir, f"{sample_name}_trimmed.fastq.gz")
        log_file = os.path.join(output_dir, f"{sample_name}_cutadapt.log")

        cutadapt_cmd = [
            "cutadapt",
            "--cores", str(threads),
            "-o", trimmed_out,
        ]
        if fwd is not None:
            cutadapt_cmd += ["-g", fwd]
        if rev is not None:
            rc_rev = _reverse_complement(rev)
            cutadapt_cmd += ["-a", rc_rev]

        cutadapt_cmd.append(nanofilt_out)

        logger.info("Running cutadapt (primer trimming) on %s", sample_name)
        logger.debug("cutadapt command: %s", " ".join(cutadapt_cmd))

        with open(log_file, "w") as log_fh:
            subprocess.run(
                cutadapt_cmd, check=True, stdout=log_fh,
                stderr=subprocess.STDOUT,
            )
        output_path = trimmed_out

    logger.info("Nanopore trimming complete. Output: %s", output_path)
    return output_path


def trim_nanopore_samples(
    samples: dict[str, str],
    output_dir: str,
    primer_set: str | None = None,
    forward_primer: str | None = None,
    reverse_primer: str | None = None,
    min_length: int = 200,
    max_length: int | None = None,
    min_quality: float = 8.0,
    threads: int = 1,
) -> dict[str, str]:
    """
    Quality-filter and optionally trim primers from multiple Nanopore samples.

    Parameters
    ----------
    samples : dict[str, str]
        Dictionary mapping sample names to single-read FASTQ paths.
    output_dir : str
        Directory where trimmed reads will be written.
    primer_set : str or None
        Built-in primer set name (see :data:`PRIMERS`).  Ignored if
        ``forward_primer`` and ``reverse_primer`` are both given.
    forward_primer : str or None
        Custom forward primer sequence.
    reverse_primer : str or None
        Custom reverse primer sequence.
    min_length : int
        Minimum read length after filtering (default: 200).
    max_length : int or None
        Maximum read length; ``None`` for no upper limit (default: None).
    min_quality : float
        Minimum mean quality score (default: 8.0).
    threads : int
        Number of threads for cutadapt (default: 1).

    Returns
    -------
    dict[str, str]
        Dictionary mapping sample names to filtered FASTQ paths.
    """
    trimmed: dict[str, str] = {}
    for sample, reads in samples.items():
        sample_out = os.path.join(output_dir, sample)
        trimmed[sample] = trim_nanopore_reads(
            reads=reads,
            output_dir=sample_out,
            min_length=min_length,
            max_length=max_length,
            min_quality=min_quality,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            primer_set=primer_set,
            threads=threads,
        )
    return trimmed


# ── Helpers ──────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgtNnYyRrSsWwKkMmBbDdHhVv",
                             "TGCAtgcaNnRrYySsWwMmKkVvHhDdBb")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (IUPAC bases supported)."""
    return seq.translate(_COMPLEMENT)[::-1]
