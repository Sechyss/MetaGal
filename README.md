# MetaGal

**MetaGal** is an open-source Python pipeline for 16S rRNA amplicon
(metabarcode) sequencing analysis.  It supports both **Illumina paired-end**
and **Oxford Nanopore MinION single-end** data, and includes a helper script
for importing data produced by the **Epi2me** software.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python ≥3.11](https://img.shields.io/badge/python-%E2%89%A53.11-blue)](https://www.python.org/)

---

## Table of contents

1. [Overview](#overview)
2. [Pipeline steps](#pipeline-steps)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick start](#quick-start)
6. [Nanopore MinION workflow](#nanopore-minion-workflow)
7. [Epi2me import](#epi2me-import)
8. [Configuration reference](#configuration-reference)
9. [Output files](#output-files)
10. [Running the tests](#running-the-tests)
11. [Citation](#citation)
12. [Licence](#licence)

---

## Overview

MetaGal wraps several established bioinformatics tools into a single,
reproducible workflow:

| Step | Tool (Illumina) | Tool (Nanopore) | Purpose |
|------|-----------------|-----------------|---------|
| 1 | FastQC / MultiQC | NanoPlot + FastQC / MultiQC | Raw read quality control |
| 2 | cutadapt | NanoFilt + cutadapt | Filtering / primer trimming |
| 3 | DADA2 paired-end (R) | DADA2 single-end (R) | Denoising & chimera removal |
| 4 | VSEARCH *or* QIIME 2 | VSEARCH *or* QIIME 2 | Taxonomic classification |
| 5 | Python (scipy / pandas) | Python (scipy / pandas) | Alpha & beta diversity |
| 6 | matplotlib | matplotlib | Figures |

All steps are configurable via a single YAML file and can be run
independently by toggling the `steps` flags in the config.

---

## Pipeline steps

```
Raw FASTQ reads
      │
      ▼
 1. Quality control
      │  Illumina: FastQC + MultiQC
      │  Nanopore: NanoPlot + FastQC + MultiQC
      ▼
 2. Trimming / Filtering
      │  Illumina: cutadapt (primer removal, paired-end)
      │  Nanopore: NanoFilt (quality + length) + optional cutadapt (primers)
      ▼
 3. Denoising (DADA2)
      │  Illumina: paired-end merge → ASV table + representative sequences
      │  Nanopore: single-end      → ASV table + representative sequences
      ▼
 4. Taxonomic classification (VSEARCH or QIIME 2)
      │  → Taxonomy table (TSV)
      ▼
 5. Diversity analysis
      │  → Alpha diversity (observed features, Shannon, Simpson)
      │  → Beta diversity  (Bray-Curtis dissimilarity matrix)
      ▼
 6. Visualisation
      │  → Alpha diversity box plot
      │  → PCoA ordination plot
      └  → Relative abundance bar chart
```

---

## Requirements

### System tools (must be on `PATH`)

| Tool | Version | Purpose | Platform |
|------|---------|---------|----------|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥0.11 | Read QC | All |
| [MultiQC](https://multiqc.info/) | ≥1.14 | QC aggregation | All |
| [cutadapt](https://cutadapt.readthedocs.io/) | ≥4.0 | Primer trimming | All |
| [R](https://www.r-project.org/) + [DADA2](https://benjjneb.github.io/dada2/) | R ≥4.2, DADA2 ≥1.26 | Denoising | All |
| [VSEARCH](https://github.com/torognes/vsearch) *(optional)* | ≥2.22 | Taxonomy (vsearch method) | All |
| [QIIME 2](https://qiime2.org/) *(optional)* | ≥2023.5 | Taxonomy (qiime2 method) | All |
| [NanoPlot](https://github.com/wdecoster/NanoPlot) | ≥1.40 | Nanopore QC | Nanopore |
| [NanoFilt](https://github.com/wdecoster/nanofilt) | ≥2.8 | Nanopore filtering | Nanopore |

### Python packages

```
cutadapt>=4.0
matplotlib>=3.7
numpy>=1.24
pandas>=2.0
PyYAML>=6.0
scikit-learn>=1.3
scipy>=1.11
```

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/Sechyss/MetaGal.git
cd MetaGal

# 2. Create and activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Install the metagal package in editable mode (optional)
pip install -e .
```

Install DADA2 in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

For Nanopore support, install NanoPlot and NanoFilt:

```bash
pip install NanoPlot NanoFilt
```

---

## Quick start

### Illumina paired-end data

1. **Copy and edit the configuration file**

   ```bash
   cp config/config.yaml my_project.yaml
   # Edit my_project.yaml – set sequencing_platform: illumina (or omit it)
   ```

2. **Prepare a sample manifest** (alternative to inline sample list)

   Create `data/manifest.tsv`:

   | sample | r1 | r2 |
   |--------|----|----|
   | sample_01 | data/raw/sample_01_R1.fastq.gz | data/raw/sample_01_R2.fastq.gz |
   | sample_02 | data/raw/sample_02_R1.fastq.gz | data/raw/sample_02_R2.fastq.gz |

3. **Run the pipeline**

   ```bash
   python pipeline.py --config my_project.yaml
   ```

---

## Nanopore MinION workflow

Use `config/config_nanopore.yaml` as your starting point:

```bash
cp config/config_nanopore.yaml my_nanopore_project.yaml
# Edit paths and parameters, then:
python pipeline.py --config my_nanopore_project.yaml
```

Key differences from the Illumina workflow:

- Set `sequencing_platform: nanopore` in the config.
- Provide a single reads file per sample using the `reads` key (instead of
  `r1` / `r2`):

  ```yaml
  samples:
    sample_01:
      reads: data/raw/sample_01.fastq.gz
  ```

  Or a TSV manifest with columns `sample` and `reads`:

  | sample | reads |
  |--------|-------|
  | sample_01 | data/raw/sample_01.fastq.gz |
  | sample_02 | data/raw/sample_02.fastq.gz |

- Trimming uses NanoFilt for quality/length filtering followed by an optional
  cutadapt primer-removal step.
- DADA2 runs in single-end mode (no paired-end merging).
- VSEARCH identity threshold defaults to `0.85` to account for higher
  Nanopore per-base error rates.

---

## Epi2me import

Oxford Nanopore's [Epi2me](https://epi2me.nanoporetech.com/) software stores
basecalled, demultiplexed reads in per-barcode directories:

```
epi2me_output/
└── fastq_pass/
    ├── barcode01/
    │   ├── reads_0.fastq.gz
    │   └── reads_1.fastq.gz
    └── barcode02/
        └── reads_0.fastq.gz
```

Use `scripts/epi2me_import.py` to concatenate the per-barcode files and
generate a manifest ready for the pipeline:

```bash
python scripts/epi2me_import.py \
    --epi2me-dir /path/to/epi2me_output \
    --sample-sheet /path/to/sample_sheet.csv \
    --output-dir data/raw \
    --manifest data/manifest_nanopore.tsv
```

The **sample sheet** is a CSV file mapping barcode IDs to sample names:

| barcode | sample_name |
|---------|-------------|
| barcode01 | patient_A |
| barcode02 | patient_B |

If `--sample-sheet` is omitted, barcode directory names are used as sample
names.  The generated `manifest_nanopore.tsv` can then be referenced in
`config_nanopore.yaml`:

```yaml
sequencing_platform: nanopore
samples: data/manifest_nanopore.tsv
```

---

## Configuration reference

See [`config/config.yaml`](config/config.yaml) for a fully annotated
Illumina example and [`config/config_nanopore.yaml`](config/config_nanopore.yaml)
for the Nanopore equivalent.  Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_dir` | – | Root directory for all outputs (**required**) |
| `samples` | – | Inline sample list or path to TSV manifest (**required**) |
| `sequencing_platform` | `illumina` | `illumina` or `nanopore` |
| `threads` | `1` | CPU threads for parallel tools |
| `metadata` | – | TSV file with sample metadata |
| `group_column` | first column | Metadata column for grouping |
| `trimming.primer_set` | `515F_806R` | Built-in primer set |
| `trimming.forward_primer` | – | Custom forward primer |
| `trimming.reverse_primer` | – | Custom reverse primer |
| `trimming.min_length` | `100` (Illumina) / `200` (Nanopore) | Minimum read length |
| `trimming.error_rate` | `0.1` | Primer mismatch rate (Illumina only) |
| `trimming.min_quality` | `8.0` | Minimum mean quality score (Nanopore only) |
| `trimming.max_length` | – | Maximum read length (Nanopore only) |
| `dada2.trunc_len_f` | `230` | Forward truncation length (Illumina) |
| `dada2.trunc_len_r` | `200` | Reverse truncation length (Illumina) |
| `dada2.trunc_len` | `0` | Truncation length, 0 = none (Nanopore) |
| `taxonomy.method` | `vsearch` | `vsearch` or `qiime2` |
| `taxonomy.database` | – | SILVA/GreenGenes FASTA or classifier .qza |
| `taxonomy.identity` | `0.97` (Illumina) / `0.85` (Nanopore) | Min alignment identity |
| `taxonomy_level` | `Phylum` | Taxonomic level for abundance plot |
| `top_n_taxa` | `10` | Top N taxa in the abundance chart |

### Built-in primer sets

| Name | Forward primer | Reverse primer | Region |
|------|----------------|----------------|--------|
| `515F_806R` | GTGYCAGCMGCCGCGGTAA | GGACTACNVGGGTWTCTAAT | V4 |
| `27F_338R` | AGAGTTTGATCMTGGCTCAG | TGCTGCCTCCCGTAGGAGT | V1-V2 |
| `341F_806R` | CCTACGGGNGGCWGCAG | GGACTACNVGGGTWTCTAAT | V3-V4 |

---

## Output files

```
results/
├── qc/
│   ├── fastqc/          # Per-sample FastQC reports
│   ├── multiqc/         # Aggregated MultiQC report
│   └── nanoplot/        # NanoPlot reports (Nanopore only)
├── trimmed/             # Filtered FASTQ files
│   └── <sample>/
├── dada2/
│   ├── asv_table.tsv    # ASV × sample count table
│   ├── rep_seqs.fasta   # Representative ASV sequences
│   └── dada2_stats.tsv  # Read tracking statistics
├── taxonomy/
│   ├── taxonomy.tsv     # Taxonomy assignments
│   └── asv_table_annotated.tsv
├── diversity/
│   ├── alpha_diversity.tsv
│   └── beta_diversity_bray_curtis.tsv
└── figures/
    ├── alpha_diversity.png
    ├── pcoa.png
    └── relative_abundance.png
```

---

## Running the tests

```bash
pip install pytest
python -m pytest tests/ -v
```

---

## Citation

If you use MetaGal in your research, please cite this repository:

> Sechyss (2024). *MetaGal: A Python pipeline for 16S metabarcode analysis*.
> GitHub. https://github.com/Sechyss/MetaGal

---

## Licence

MetaGal is released under the **GNU General Public License v3.0**.
See [LICENSE](LICENSE) for the full text.