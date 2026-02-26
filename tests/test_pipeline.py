"""
Tests for pipeline.py – configuration loading and sample parsing.
"""

import textwrap

import pandas as pd
import pytest

from pipeline import load_config, load_samples


# ── load_config ───────────────────────────────────────────────────────────────

def test_load_config_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(str(tmp_path / "nonexistent.yaml"))


def test_load_config_missing_keys(tmp_path):
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("threads: 4\n")
    with pytest.raises(ValueError, match="Missing required configuration keys"):
        load_config(str(cfg))


def test_load_config_valid(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        textwrap.dedent("""\
            output_dir: /tmp/results
            samples:
              s1:
                r1: /data/s1_R1.fastq.gz
                r2: /data/s1_R2.fastq.gz
        """)
    )
    config = load_config(str(cfg))
    assert config["output_dir"] == "/tmp/results"
    assert "s1" in config["samples"]


# ── load_samples ──────────────────────────────────────────────────────────────

def test_load_samples_from_dict():
    samples_cfg = {
        "sample_A": {"r1": "A_R1.fastq.gz", "r2": "A_R2.fastq.gz"},
        "sample_B": {"r1": "B_R1.fastq.gz", "r2": "B_R2.fastq.gz"},
    }
    result = load_samples(samples_cfg)
    assert set(result.keys()) == {"sample_A", "sample_B"}
    assert result["sample_A"] == ("A_R1.fastq.gz", "A_R2.fastq.gz")


def test_load_samples_from_tsv(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample\tr1\tr2\nS1\tS1_R1.fq\tS1_R2.fq\n")
    result = load_samples(str(manifest))
    assert "S1" in result
    assert result["S1"] == ("S1_R1.fq", "S1_R2.fq")
