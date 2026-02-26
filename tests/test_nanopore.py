"""
Tests for Nanopore support – pipeline sample loading, trimming helpers,
and the Epi2me import script.
"""

from __future__ import annotations

import gzip
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from pipeline import load_samples_nanopore
from metagal.trimming import trim_nanopore_samples, PRIMERS
from scripts.epi2me_import import (
    find_barcode_dirs,
    collect_fastq_files,
    concatenate_fastq,
    load_sample_sheet,
    import_epi2me,
)


# ── load_samples_nanopore ─────────────────────────────────────────────────────

def test_load_samples_nanopore_from_dict():
    cfg = {
        "sample_A": {"reads": "A.fastq.gz"},
        "sample_B": {"reads": "B.fastq.gz"},
    }
    result = load_samples_nanopore(cfg)
    assert set(result.keys()) == {"sample_A", "sample_B"}
    assert result["sample_A"] == "A.fastq.gz"
    assert result["sample_B"] == "B.fastq.gz"


def test_load_samples_nanopore_from_tsv(tmp_path):
    manifest = tmp_path / "manifest_nanopore.tsv"
    manifest.write_text("sample\treads\nS1\tS1.fastq.gz\nS2\tS2.fastq.gz\n")
    result = load_samples_nanopore(str(manifest))
    assert set(result.keys()) == {"S1", "S2"}
    assert result["S1"] == "S1.fastq.gz"


# ── trim_nanopore_samples (ValueError for unknown primer set) ─────────────────

def test_trim_nanopore_unknown_primer_set(tmp_path):
    """Requesting an unrecognised primer set should raise ValueError."""
    # Create a dummy reads file so the FileNotFoundError is not raised first
    reads = tmp_path / "sample.fastq.gz"
    reads.write_bytes(b"")

    from metagal.trimming import trim_nanopore_reads
    with pytest.raises(ValueError):
        trim_nanopore_reads(
            reads=str(reads),
            output_dir=str(tmp_path / "out"),
            primer_set="UNKNOWN_SET",
        )


# ── find_barcode_dirs ─────────────────────────────────────────────────────────

def test_find_barcode_dirs_top_level(tmp_path):
    (tmp_path / "barcode01").mkdir()
    (tmp_path / "barcode02").mkdir()
    (tmp_path / "not_a_barcode").mkdir()

    result = find_barcode_dirs(str(tmp_path))
    assert "barcode01" in result
    assert "barcode02" in result
    assert "not_a_barcode" not in result


def test_find_barcode_dirs_fastq_pass(tmp_path):
    fq_pass = tmp_path / "fastq_pass"
    fq_pass.mkdir()
    (fq_pass / "barcode03").mkdir()
    (fq_pass / "barcode04").mkdir()

    result = find_barcode_dirs(str(tmp_path))
    assert "barcode03" in result
    assert "barcode04" in result


def test_find_barcode_dirs_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_barcode_dirs(str(tmp_path / "nonexistent"))


def test_find_barcode_dirs_case_insensitive(tmp_path):
    (tmp_path / "Barcode01").mkdir()
    result = find_barcode_dirs(str(tmp_path))
    assert "barcode01" in result


# ── collect_fastq_files ───────────────────────────────────────────────────────

def test_collect_fastq_files(tmp_path):
    (tmp_path / "reads.fastq").write_text("@r1\nACGT\n+\nIIII\n")
    (tmp_path / "reads2.fastq.gz").write_bytes(b"")
    (tmp_path / "other.txt").write_text("ignore me")

    files = collect_fastq_files(tmp_path)
    names = [f.name for f in files]
    assert "reads.fastq" in names
    assert "reads2.fastq.gz" in names
    assert "other.txt" not in names


# ── concatenate_fastq ─────────────────────────────────────────────────────────

def test_concatenate_fastq_plain(tmp_path):
    fq = tmp_path / "reads.fastq"
    content = b"@r1\nACGT\n+\nIIII\n"
    fq.write_bytes(content)

    out = tmp_path / "merged.fastq.gz"
    concatenate_fastq([fq], out)

    with gzip.open(out, "rb") as fh:
        assert fh.read() == content


def test_concatenate_fastq_gzipped(tmp_path):
    fq = tmp_path / "reads.fastq.gz"
    content = b"@r1\nACGT\n+\nIIII\n"
    with gzip.open(fq, "wb") as fh:
        fh.write(content)

    out = tmp_path / "merged.fastq.gz"
    concatenate_fastq([fq], out)

    with gzip.open(out, "rb") as fh:
        assert fh.read() == content


def test_concatenate_fastq_multiple(tmp_path):
    fq1 = tmp_path / "a.fastq"
    fq2 = tmp_path / "b.fastq"
    fq1.write_bytes(b"@r1\nAAAA\n+\nIIII\n")
    fq2.write_bytes(b"@r2\nCCCC\n+\nIIII\n")

    out = tmp_path / "merged.fastq.gz"
    concatenate_fastq([fq1, fq2], out)

    with gzip.open(out, "rb") as fh:
        merged = fh.read()
    assert b"AAAA" in merged
    assert b"CCCC" in merged


def test_concatenate_fastq_empty_list(tmp_path):
    with pytest.raises(ValueError):
        concatenate_fastq([], tmp_path / "out.fastq.gz")


# ── load_sample_sheet ─────────────────────────────────────────────────────────

def test_load_sample_sheet_basic(tmp_path):
    csv = tmp_path / "sheet.csv"
    csv.write_text("barcode,sample_name\nbarcode01,patient_A\nbarcode02,patient_B\n")
    mapping = load_sample_sheet(str(csv))
    assert mapping["barcode01"] == "patient_A"
    assert mapping["barcode02"] == "patient_B"


def test_load_sample_sheet_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_sample_sheet(str(tmp_path / "nonexistent.csv"))


def test_load_sample_sheet_missing_column(tmp_path):
    csv = tmp_path / "sheet.csv"
    csv.write_text("wrong_col,sample_name\nbarcode01,patient_A\n")
    with pytest.raises(KeyError):
        load_sample_sheet(str(csv))


def test_load_sample_sheet_custom_cols(tmp_path):
    csv = tmp_path / "sheet.csv"
    csv.write_text("bc,name\nbarcode01,S1\n")
    mapping = load_sample_sheet(str(csv), barcode_col="bc", name_col="name")
    assert mapping["barcode01"] == "S1"


# ── import_epi2me ─────────────────────────────────────────────────────────────

def _make_epi2me_dir(tmp_path: Path) -> tuple[Path, Path]:
    """Create a minimal Epi2me-like output directory structure."""
    fq_pass = tmp_path / "fastq_pass"
    bc01 = fq_pass / "barcode01"
    bc02 = fq_pass / "barcode02"
    bc01.mkdir(parents=True)
    bc02.mkdir(parents=True)

    (bc01 / "reads_0.fastq").write_bytes(b"@r1\nACGT\n+\nIIII\n")
    (bc02 / "reads_0.fastq").write_bytes(b"@r2\nTTTT\n+\nIIII\n")

    return fq_pass, tmp_path


def test_import_epi2me_no_sample_sheet(tmp_path):
    _, epi2me_root = _make_epi2me_dir(tmp_path)
    out_dir = tmp_path / "raw"
    manifest = tmp_path / "manifest.tsv"

    df = import_epi2me(
        epi2me_dir=str(epi2me_root),
        output_dir=str(out_dir),
        manifest_path=str(manifest),
    )

    assert len(df) == 2
    assert set(df["sample"].tolist()) == {"barcode01", "barcode02"}
    assert manifest.exists()
    # Each output FASTQ should exist
    for reads_path in df["reads"]:
        assert Path(reads_path).exists()


def test_import_epi2me_with_sample_sheet(tmp_path):
    _, epi2me_root = _make_epi2me_dir(tmp_path)
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("barcode,sample_name\nbarcode01,patient_A\nbarcode02,patient_B\n")
    out_dir = tmp_path / "raw"
    manifest = tmp_path / "manifest.tsv"

    df = import_epi2me(
        epi2me_dir=str(epi2me_root),
        output_dir=str(out_dir),
        manifest_path=str(manifest),
        sample_sheet_path=str(sheet),
    )

    assert set(df["sample"].tolist()) == {"patient_A", "patient_B"}


def test_import_epi2me_manifest_columns(tmp_path):
    _, epi2me_root = _make_epi2me_dir(tmp_path)
    out_dir = tmp_path / "raw"
    manifest = tmp_path / "manifest.tsv"

    import_epi2me(
        epi2me_dir=str(epi2me_root),
        output_dir=str(out_dir),
        manifest_path=str(manifest),
    )

    loaded = pd.read_csv(manifest, sep="\t")
    assert "sample" in loaded.columns
    assert "reads" in loaded.columns


def test_import_epi2me_missing_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        import_epi2me(
            epi2me_dir=str(tmp_path / "nonexistent"),
            output_dir=str(tmp_path / "raw"),
            manifest_path=str(tmp_path / "manifest.tsv"),
        )


def test_import_epi2me_no_barcodes(tmp_path):
    (tmp_path / "some_other_dir").mkdir()
    with pytest.raises(RuntimeError, match="No barcode directories"):
        import_epi2me(
            epi2me_dir=str(tmp_path),
            output_dir=str(tmp_path / "raw"),
            manifest_path=str(tmp_path / "manifest.tsv"),
        )
