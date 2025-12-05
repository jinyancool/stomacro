# Integrated Spatial Multi-Omics Analysis of Resectable Oral Squamous Cell Carcinoma (OSCC)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/) [![Nextflow](https://img.shields.io/badge/nextflow-%3E%2023.10-099acd.svg)](https://www.nextflow.io) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32-brightgreen.svg)](https://snakemake.readthedocs.io)

**Repository**: https://github.com/jinyancool/stomacro

**Principal Investigator**: Jhuanglab

**Contact**: hiekeen \[at\] gmail.com

## Project Overview

**StomaCro** (Spatial Tumor Microenvironment analysis of Oral Cancer) is a comprehensive, reproducible research compendium for the study:

> **"Integrated spatial multi-omics defines spatially resolved molecular subtypes and actionable therapeutic vulnerabilities in resectable oral squamous cell carcinoma"**

Oral squamous cell carcinoma (OSCC) is the most common head-and-neck malignancy in East Asia, with \>30–40% recurrence rate after curative surgery. Bulk and dissociated single-cell omics ignore the critical spatial organization of the tumor microenvironment (TME). This project closes that gap by integrating four state-of-the-art **spatial multi-omics platforms** on the same set of resectable (stage I–III) OSCC samples from an East Asian cohort:

-   10x Genomics Visium (whole-transcriptome spatial transcriptomics)
-   10x Genomics Xenium In Situ (subcellular-resolution targeted RNA panels)
-   NanoString CosMx Spatial Molecular Imager (1,000-plex RNA + 64-plex protein, single-cell)
-   Akoya CODEX (50+-plex highly validated protein panel)

### Key Scientific Aims

1.  Define spatially restricted molecular subtypes of OSCC
2.  Dissect tumor core vs. invasive front heterogeneity (EMT, hypoxia, immune exclusion)
3.  Map cell–cell interaction networks (e.g., CAF → tumor polyamine hubs, CXCL12–CXCR4 axes)
4.  Identify spatially confined therapeutic targets (PARP inhibitors, CXCR4 blockade, bispecific antibodies)
5.  Build AI-enhanced prognostic models incorporating spatial features

## Quick Start 

``` bash
# 1. Clone the repository
git clone https://github.com/jinyancool/stomacro.git
cd stomacro
```

## Software Environment

-   Python ≥ 3.9
-   R ≥ 4.3
-   Scanpy, Seurat, Squidpy, Giotto, STutility
-   CellPhoneDB v4, LIANA, NicheNet (spatial mode)
-   TensorFlow / PyTorch for deep-learning spatial modeling

Detailed environments are in `envs/`.

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Funded by \[Your Funding Agency/Grant Number\]. We thank 10x Genomics, NanoString, and Akoya Biosciences for technical support and early-access programs.
