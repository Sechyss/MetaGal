"""
Denoising module for 16S metabarcode analysis.

Wraps DADA2 (via an R subprocess) to learn error rates, denoise reads,
merge paired-end reads (Illumina), or process single-end long reads
(Nanopore), and remove chimeras, producing an ASV table.
"""

import json
import logging
import os
import subprocess
import textwrap
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Inline R script executed by run_dada2()
_DADA2_SCRIPT = textwrap.dedent("""\
    library(dada2)
    library(jsonlite)
    library(Biostrings)

    args        <- commandArgs(trailingOnly = TRUE)
    input_json  <- args[1]
    output_dir  <- args[2]
    trunc_len_f <- as.integer(args[3])
    trunc_len_r <- as.integer(args[4])
    threads     <- as.integer(args[5])

    params <- fromJSON(input_json)   # list: sample -> c(R1, R2)
    fwd_files <- sapply(params, `[[`, 1)
    rev_files <- sapply(params, `[[`, 2)
    sample_names <- names(params)

    # Learn error rates
    err_fwd <- learnErrors(fwd_files, multithread = threads)
    err_rev <- learnErrors(rev_files, multithread = threads)

    # Dereplicate
    derep_fwd <- derepFastq(fwd_files)
    derep_rev <- derepFastq(rev_files)
    names(derep_fwd) <- sample_names
    names(derep_rev) <- sample_names

    # Denoise
    dada_fwd <- dada(derep_fwd, err = err_fwd, multithread = threads)
    dada_rev <- dada(derep_rev, err = err_rev, multithread = threads)

    # Merge paired reads
    merged <- mergePairs(dada_fwd, derep_fwd, dada_rev, derep_rev)

    # Build sequence table
    seqtab <- makeSequenceTable(merged)

    # Remove chimeras
    seqtab_nochim <- removeBimeraDenovo(
        seqtab, method = "consensus", multithread = threads
    )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Write ASV table (samples as rows, ASVs as columns)
    asv_table_path <- file.path(output_dir, "asv_table.tsv")
    write.table(
        t(seqtab_nochim),
        file = asv_table_path,
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    # Write representative sequences (FASTA)
    fasta_path <- file.path(output_dir, "rep_seqs.fasta")
    seqs <- colnames(seqtab_nochim)
    asv_ids <- paste0("ASV", seq_along(seqs))
    writeXStringSet(
        DNAStringSet(setNames(seqs, asv_ids)),
        filepath = fasta_path
    )

    # Write summary stats
    stats <- data.frame(
        sample   = sample_names,
        input    = sapply(dada_fwd, function(x) sum(getUniques(x))),
        merged   = sapply(merged,   function(x) sum(getUniques(x))),
        nonchim  = rowSums(seqtab_nochim)
    )
    write.table(stats,
        file = file.path(output_dir, "dada2_stats.tsv"),
        sep = "\\t", quote = FALSE, row.names = FALSE
    )

    message("DADA2 complete. ASV table: ", asv_table_path)
""")


def run_dada2(
    samples: dict[str, tuple[str, str]],
    output_dir: str,
    trunc_len_f: int = 230,
    trunc_len_r: int = 200,
    threads: int = 1,
) -> tuple[str, str]:
    """
    Denoise paired-end reads with DADA2 and remove chimeras.

    Parameters
    ----------
    samples : dict[str, tuple[str, str]]
        Dictionary mapping sample names to (R1, R2) trimmed FASTQ paths.
    output_dir : str
        Directory where DADA2 outputs will be written.
    trunc_len_f : int
        Truncate forward reads at this position (default: 230).
    trunc_len_r : int
        Truncate reverse reads at this position (default: 200).
    threads : int
        Number of threads for DADA2 (default: 1).

    Returns
    -------
    tuple[str, str]
        Paths to the ASV table TSV and representative sequences FASTA.

    Raises
    ------
    FileNotFoundError
        If any input FASTQ file does not exist.
    subprocess.CalledProcessError
        If the R script exits with a non-zero return code.
    """
    for sample, (r1, r2) in samples.items():
        for f in (r1, r2):
            if not Path(f).exists():
                raise FileNotFoundError(f"Input file not found: {f}")

    os.makedirs(output_dir, exist_ok=True)

    # Write R script to a temp file
    script_path = os.path.join(output_dir, "_dada2_run.R")
    with open(script_path, "w") as fh:
        fh.write(_DADA2_SCRIPT)

    # Serialise sample paths to JSON for the R script
    input_json = json.dumps({s: list(paths) for s, paths in samples.items()})

    cmd = [
        "Rscript", script_path,
        input_json,
        output_dir,
        str(trunc_len_f),
        str(trunc_len_r),
        str(threads),
    ]

    logger.info("Running DADA2 on %d sample(s)", len(samples))
    logger.debug("Command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    asv_table = os.path.join(output_dir, "asv_table.tsv")
    rep_seqs = os.path.join(output_dir, "rep_seqs.fasta")
    logger.info("DADA2 complete. ASV table: %s", asv_table)
    return asv_table, rep_seqs


# Inline R script for single-end Nanopore reads (no paired merging)
_DADA2_NANOPORE_SCRIPT = textwrap.dedent("""\
    library(dada2)
    library(jsonlite)
    library(Biostrings)

    args        <- commandArgs(trailingOnly = TRUE)
    input_json  <- args[1]
    output_dir  <- args[2]
    trunc_len   <- as.integer(args[3])   # 0 = no truncation
    threads     <- as.integer(args[4])

    params <- fromJSON(input_json)   # list: sample -> reads path
    read_files   <- unlist(params)
    sample_names <- names(params)

    # Learn error rates from single-end reads
    err <- learnErrors(read_files, multithread = threads)

    # Dereplicate
    derep <- derepFastq(read_files)
    names(derep) <- sample_names

    # Denoise (single-end – no merging step)
    dada_out <- dada(derep, err = err, multithread = threads)

    # Build sequence table
    if (trunc_len > 0) {
        seqtab <- makeSequenceTable(dada_out)
        seqtab <- seqtab[, nchar(colnames(seqtab)) == trunc_len]
    } else {
        seqtab <- makeSequenceTable(dada_out)
    }

    # Remove chimeras
    seqtab_nochim <- removeBimeraDenovo(
        seqtab, method = "consensus", multithread = threads
    )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Write ASV table (samples as rows, ASVs as columns)
    asv_table_path <- file.path(output_dir, "asv_table.tsv")
    write.table(
        t(seqtab_nochim),
        file = asv_table_path,
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    # Write representative sequences (FASTA)
    fasta_path <- file.path(output_dir, "rep_seqs.fasta")
    seqs    <- colnames(seqtab_nochim)
    asv_ids <- paste0("ASV", seq_along(seqs))
    writeXStringSet(
        DNAStringSet(setNames(seqs, asv_ids)),
        filepath = fasta_path
    )

    # Write summary stats
    stats <- data.frame(
        sample  = sample_names,
        input   = sapply(dada_out, function(x) sum(getUniques(x))),
        nonchim = rowSums(seqtab_nochim)
    )
    write.table(stats,
        file = file.path(output_dir, "dada2_stats.tsv"),
        sep = "\\t", quote = FALSE, row.names = FALSE
    )

    message("DADA2 (Nanopore) complete. ASV table: ", asv_table_path)
""")


def run_dada2_nanopore(
    samples: dict[str, str],
    output_dir: str,
    trunc_len: int = 0,
    threads: int = 1,
) -> tuple[str, str]:
    """
    Denoise single-end Nanopore 16S reads with DADA2 and remove chimeras.

    This function runs DADA2 in single-end mode (no paired-end merging),
    which is appropriate for long-read Nanopore amplicon data.

    Parameters
    ----------
    samples : dict[str, str]
        Dictionary mapping sample names to single-read (FASTQ) paths.
    output_dir : str
        Directory where DADA2 outputs will be written.
    trunc_len : int
        Truncate reads to exactly this length; ``0`` disables truncation
        (recommended for variable-length Nanopore reads, default: 0).
    threads : int
        Number of threads for DADA2 (default: 1).

    Returns
    -------
    tuple[str, str]
        Paths to the ASV table TSV and representative sequences FASTA.

    Raises
    ------
    FileNotFoundError
        If any input FASTQ file does not exist.
    subprocess.CalledProcessError
        If the R script exits with a non-zero return code.
    """
    for sample, reads in samples.items():
        if not Path(reads).exists():
            raise FileNotFoundError(f"Input file not found: {reads}")

    os.makedirs(output_dir, exist_ok=True)

    script_path = os.path.join(output_dir, "_dada2_nanopore_run.R")
    with open(script_path, "w") as fh:
        fh.write(_DADA2_NANOPORE_SCRIPT)

    input_json = json.dumps({s: reads for s, reads in samples.items()})

    cmd = [
        "Rscript", script_path,
        input_json,
        output_dir,
        str(trunc_len),
        str(threads),
    ]

    logger.info("Running DADA2 (Nanopore) on %d sample(s)", len(samples))
    logger.debug("Command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    asv_table = os.path.join(output_dir, "asv_table.tsv")
    rep_seqs = os.path.join(output_dir, "rep_seqs.fasta")
    logger.info("DADA2 (Nanopore) complete. ASV table: %s", asv_table)
    return asv_table, rep_seqs


def load_asv_table(asv_table_path: str) -> pd.DataFrame:
    """
    Load an ASV table produced by :func:`run_dada2`.

    Parameters
    ----------
    asv_table_path : str
        Path to the TSV file produced by DADA2.

    Returns
    -------
    pd.DataFrame
        DataFrame with ASVs as rows and samples as columns.

    Raises
    ------
    FileNotFoundError
        If the ASV table file does not exist.
    """
    if not Path(asv_table_path).exists():
        raise FileNotFoundError(f"ASV table not found: {asv_table_path}")

    return pd.read_csv(asv_table_path, sep="\t", index_col=0)
