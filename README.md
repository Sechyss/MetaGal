# MetaGal

**MetaGal** is an open-source Python pipeline for 16S rRNA amplicon
(metabarcode) sequencing analysis.  It takes paired-end Illumina FASTQ
files from raw reads through to alpha/beta diversity metrics and
publication-ready figures.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python ≥3.11](https://img.shields.io/badge/python-%E2%89%A53.11-blue)](https://www.python.org/)

---

## Table of contents

1. [Overview](#overview)
2. [Pipeline steps](#pipeline-steps)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick start](#quick-start)
6. [Configuration reference](#configuration-reference)
7. [Output files](#output-files)
8. [Running the tests](#running-the-tests)
9. [Citation](#citation)
10. [Licence](#licence)

---

## Overview

MetaGal wraps several established bioinformatics tools into a single,
reproducible workflow:

| Step | Tool | Purpose |
|------|------|---------|
| 1 | FastQC / MultiQC | Raw read quality control |
| 2 | cutadapt | Primer trimming |
| 3 | DADA2 (R) | Denoising, paired-end merging, chimera removal |
| 4 | VSEARCH *or* QIIME 2 | Taxonomic classification |
| 5 | Python (scipy / pandas) | Alpha & beta diversity |
| 6 | matplotlib | Figures |

All steps are configurable via a single YAML file and can be run
independently by toggling the `steps` flags in the config.

---

## Pipeline steps

```
Raw FASTQ reads
      │
      ▼
 1. Quality control (FastQC + MultiQC)
      │
      ▼
 2. Primer trimming (cutadapt)
      │
      ▼
 3. Denoising (DADA2)
      │  → ASV table (TSV)
      │  → Representative sequences (FASTA)
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

| Tool | Version | Purpose |
|------|---------|---------|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥0.11 | Read QC |
| [MultiQC](https://multiqc.info/) | ≥1.14 | QC aggregation |
| [cutadapt](https://cutadapt.readthedocs.io/) | ≥4.0 | Primer trimming |
| [R](https://www.r-project.org/) + [DADA2](https://benjjneb.github.io/dada2/) | R ≥4.2, DADA2 ≥1.26 | Denoising |
| [VSEARCH](https://github.com/torognes/vsearch) *(optional)* | ≥2.22 | Taxonomy (vsearch method) |
| [QIIME 2](https://qiime2.org/) *(optional)* | ≥2023.5 | Taxonomy (qiime2 method) |

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

---

## Quick start

1. **Copy and edit the configuration file**

   ```bash
   cp config/config.yaml my_project.yaml
   # Edit my_project.yaml to point to your FASTQ files and database
   ```

2. **Prepare a sample manifest** (alternative to inline sample list)

   Create `data/manifest.tsv` with the following columns:

   | sample | r1 | r2 |
   |--------|----|----|
   | sample_01 | data/raw/sample_01_R1.fastq.gz | data/raw/sample_01_R2.fastq.gz |
   | sample_02 | data/raw/sample_02_R1.fastq.gz | data/raw/sample_02_R2.fastq.gz |

3. **Prepare a metadata file** (`data/metadata.tsv`)

   | sample_id | treatment | timepoint |
   |-----------|-----------|-----------|
   | sample_01 | control | 0 |
   | sample_02 | treated | 24 |

4. **Run the pipeline**

   ```bash
   python pipeline.py --config my_project.yaml
   ```

   Run a specific step only (e.g. diversity + visualisation):

   ```yaml
   # In my_project.yaml, set unwanted steps to false:
   steps:
     quality_control: false
     trimming: false
     denoising: false
     taxonomy: false
     diversity: true
     visualisation: true

   # Provide pre-computed ASV table when skipping denoising
   asv_table: results/dada2/asv_table.tsv
   rep_seqs:  results/dada2/rep_seqs.fasta
   ```

---

## Configuration reference

See [`config/config.yaml`](config/config.yaml) for a fully annotated
example.  Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_dir` | – | Root directory for all outputs (**required**) |
| `samples` | – | Inline sample list or path to TSV manifest (**required**) |
| `threads` | `1` | CPU threads for parallel tools |
| `metadata` | – | TSV file with sample metadata |
| `group_column` | first column | Metadata column for grouping |
| `trimming.primer_set` | `515F_806R` | Built-in primer set (V4 region) |
| `trimming.forward_primer` | – | Custom forward primer (overrides `primer_set`) |
| `trimming.reverse_primer` | – | Custom reverse primer |
| `trimming.min_length` | `100` | Minimum read length after trimming |
| `dada2.trunc_len_f` | `230` | Forward read truncation length |
| `dada2.trunc_len_r` | `200` | Reverse read truncation length |
| `taxonomy.method` | `vsearch` | `vsearch` or `qiime2` |
| `taxonomy.database` | – | SILVA/GreenGenes FASTA (vsearch) or classifier .qza (qiime2) |
| `taxonomy.identity` | `0.97` | Minimum alignment identity (vsearch) |
| `taxonomy_level` | `Phylum` | Taxonomic level for relative abundance plot |
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
│   └── multiqc/         # Aggregated MultiQC report
├── trimmed/             # Primer-trimmed FASTQ files
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