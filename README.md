# An atlas of colonization factors in the human gut microbiome reveals ecological strategies and inflammatory bowel disease signatures

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21643081-blue)](https://doi.org/10.5281/zenodo.21643081)
[![R](https://img.shields.io/badge/R-4.4.3-276DC3)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.10.16-3776AB)](https://www.python.org/)

This repository contains the analysis code, processed intermediate data, randomization analyses, and figure-generation scripts associated with the study:

> **Meng, Jin-Xin, et al.** *An atlas of colonization factors in the human gut microbiome reveals ecological strategies and inflammatory bowel disease signatures.* Nature Communications (2026).

The project systematically maps gut microbial colonization factors (CFs) across the Unified Human Gastrointestinal Protein/Genome resources and evaluates their ecological organization, taxonomic distribution, disease-associated shifts, predictive potential, and metatranscriptomic activity in inflammatory bowel disease (IBD).

---

## Overview

The analysis is organized around four major components:

1. **Construction of a gut microbial CF atlas**
   - Identification of CF homologs in the UHGP protein catalogue.
   - Propagation of CF annotations from UHGP-90 representatives to UHGG genomes/species.
   - Taxonomic and phylogenetic characterization of CF repertoires.

2. **Ecological organization of CF repertoires**
   - Presence/absence profiling across gut microbial species.
   - Jaccard-distance analyses and phylogenetic/taxonomic structure.
   - PAM clustering of species according to CF repertoires.
   - Functional interpretation using KEGG Orthology profiles.

3. **IBD metagenomic analysis**
   - Quantification of CF genes in multiple public shotgun metagenomic cohorts.
   - Species-level community profiling using the UHGG reference database.
   - Diversity, compositional, concordance, differential-abundance, meta-analysis, and classification analyses.

4. **Metatranscriptomic analysis**
   - Quantification of CF transcription in IBD metatranscriptomes.
   - Comparison of CF abundance and expression.
   - Disease-associated transcriptional signatures and predictive analyses.

---

## Repository structure

```text
25-CF-IBD-MGX/
├── README.md
├── pipeline/                # CF identification, MGX/MTX profiling and reference files
├── data/                    # Processed cohort-level profiles and sample metadata
├── randomization/           # Randomization controls, scripts and intermediate results
├── scripts/                 # Shared R/Python helper functions
├── Figure1/                 # Scripts and source data for Figure 1
├── Figure2/                 # Scripts and source data for Figure 2
├── Figure3/                 # Scripts and source data for Figure 3
├── Figure4/                 # Scripts and source data for Figure 4
├── Figure5/                 # Scripts and source data for Figure 5
├── Figure6/                 # Scripts and source data for Figure 6
└── Figure7/                 # Scripts and source data for Figure 7
```

Most large tabular files are compressed with `bzip2` (`*.bz2`).

---

## Public cohorts

| Cohort | BioProject |
|---|---|
| BushmanFD_2020 | PRJNA562600 |
| FranzosaEA_2018 | PRJNA400072 |
| HallAB_2017 | PRJNA385949 |
| HeQ_2017 | PRJEB15371 |
| KumbhariA_2024 | PRJNA993675 |
| LloydPriceJ_2019 | PRJNA398089 |
| SchirmerM_2018 | PRJNA389280 |
| SchirmerM_2024 | PRJNA436359 |
| WengY_2019 | PRJNA429990 |
| YanQ_2023c | PRJEB67456 |

Metatranscriptomic analyses use: LloydPriceJ_2019 and SchirmerM_2018.

---

## Shared analysis functions

| Script | Main purpose |
|---|---|
| `calcu_difference.R` | Group-wise statistical comparisons |
| `calcu_metafor.R` | Meta-analysis utilities |
| `calcu_metafor-0.1.1.R` | Extended/versioned meta-analysis functions |
| `diversity.R` | Diversity-related calculations |
| `model_randomforest.R` | Random-forest modelling and cross-validation |
| `palette.R` | Shared color palettes |
| `plot_Procrustes.R` | Procrustes analysis/visualization |
| `plot_pie.R` | Pie/hemisphere-pie plotting |
| `plot_roc.R` | ROC visualization |
| `run-fast-combine-matrix.py` | Combine per-sample output files into matrices |
| `transform_rc.R` | Read-count/abundance profile transformation |

---

## Software environment

Analyses were performed under:

- Ubuntu 24.04.4
- Python 3.10.16
- R 4.4.3

### Command-line software

| Software | Version / notes |
|---|---|
| DIAMOND | 2.1.11.165 |
| eggNOG-mapper | 2.1 |
| eggNOG database | 5.0 |
| PhyloPhlAn | 3.1.68 |
| MCL | 22-282 |
| fastp | 0.24 |
| Bowtie2 | 2.5.4 |
| Samtools | 1.21 |
| Kraken2 | 2.1.3 |
| Bracken | 3.0.1 |
| Cytoscape | 3.10.3 |
| iTOL | 7.3 |

### Major R packages

| Package | Version |
|---|---|
| clusterSim | 0.51-5 |
| vegan | 2.7 |
| ape | 5.8 |
| metafor | 3.8-0 |
| randomForest | 4.7 |
| pROC | 1.19 |
| rstatix | 0.7.2 |
| DESeq2 | 1.46 |
| clusterProfiler | 4.17 |
| ggplot2 | 3.5.2 |
| ggtern | 3.5 |
| ggpubr | 0.6.1 |

### Major Python packages

| Package | Version |
|---|---|
| NumPy | 2.1.3 |
| pandas | 2.2.3 |

---

## External resources

This project relies on external reference resources that are not fully redistributed here:

- genotype-habitat association (GHA) project / CF reference proteins;
- UHGP / UHGP-90 protein catalogues;
- UHGG v2.0.2 genome catalogue and metadata;
- UHGG Kraken2 database;
- eggNOG v5.0;
- public sequencing datasets listed above.

---

## Citation

If you use this repository, CF catalogue, code, or processed data, please cite:

> Meng, Jin-Xin, et al. **An atlas of colonization factors in the human gut microbiome reveals ecological strategies and inflammatory bowel disease signatures.** *Nature Communications* (2026).

Archived code/data release:

> Zenodo. https://doi.org/10.5281/zenodo.21643081

---

## Acknowledgements

This repository incorporates analyses based on public metagenomic and metatranscriptomic datasets, external reference resources, and open-source bioinformatics software. We thank the authors and developers of the genotype-habitat association (GHA) project, UHGP/UHGG, and the other tools and resources used in this study. Please cite the corresponding original resources, software, and cohort publications when reusing these data or workflows.
