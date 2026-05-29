# WORMS

This is part of my PhD project, which involved originally 2 helminth species for the study of comple life cycles.
*Camallanus lacustris* and *Schistocephalus solidus*.
In the Copepod folder, analyses related to the first host of both.

In the Schistocephalus folder: analysis pipelines, scripts, and resources for studying gene expression and life-stage decoupling in the tapeworm *Schistocephalus solidus*.

Repository: [LauraGramolini/Worms – Schistocephalus directory](https://github.com/LauraGramolini/Worms/tree/master/Schistocephalus?utm_source=chatgpt.com)

---

## Overview

This repository contains scripts and analysis workflows used for transcriptomic analyses of the cestode parasite *Schistocephalus solidus* across its complex life cycle. The project focuses on understanding how gene expression changes across hosts and functional stages, with particular emphasis on:

* Adaptive decoupling across life stages
* Host-specific versus stage-specific gene expression
* Functional enrichment and GO term analyses
* Differential expression across developmental transitions
* Dual RNA-seq host–parasite analyses

The analyses support research on phenotypic plasticity, host transitions, and evolutionary adaptation in parasites with complex life cycles.

---

## Biological System

*Schistocephalus solidus* is a tapeworm with a complex life cycle involving multiple hosts:

1. **Copepod** (first intermediate host)
2. **Stickleback fish** (second intermediate host)
3. **Bird** (definitive host)

The repository includes analyses spanning several functional stages including:

* Free-living larval stages
* Infection stages
* Growth stages
* Transmission stages
* Reproductive stages

More information about the organism can be found at:

* [WormBase ParaSite – Schistocephalus solidus](https://parasite.wormbase.org/Schistocephalus_solidus_prjeb527/Info/Index/?utm_source=chatgpt.com)
* Schistocephalus solidus

---

## Repository Structure

Example structure (may vary slightly depending on updates):

```text

Schistocephalus/
│
├── Bash_scripts/
├── R/
├── data/
└── README.md

```

---

## Main Analyses

### Differential Gene Expression

The repository includes workflows for:

* Read preprocessing
* Transcript quantification
* Differential expression analysis
* Stage-to-stage comparisons
* Host-transition comparisons

Typical tools used in the pipeline may include:

* STAR
* featureCounts
* DESeq2
* Trinity
* eggNOG-mapper

---

### Functional Enrichment

Gene Ontology (GO) enrichment analyses are used to identify:

* Stage-specific biological functions
* Shared versus decoupled pathways
* Concordant and discordant enrichment patterns

The analyses investigate whether biologically similar stages reuse the same genes or rely on independent regulatory programs.

---

### Correlation and Decoupling Analyses

A major focus of the project is testing the adaptive decoupling hypothesis through:

* Gene-level expression correlations
* Shared DEG analyses
* Functional overlap analyses
* Consecutive versus non-consecutive stage comparisons

These analyses evaluate how independently different life stages evolve and regulate transcription.

---

## Requirements

The repository primarily uses:

* R
* Bash
* Python

Suggested R packages:

```r
DESeq2
edgeR
ggplot2
clusterProfiler
GOSemSim
tidyverse
```

Suggested Python packages:

```python
pandas
numpy
matplotlib
seaborn
scipy
```

---

## Usage

### Clone the repository

```bash
git clone https://github.com/LauraGramolini/Worms.git
cd Worms/Schistocephalus
```

### Run analyses

Example workflow:

```bash
Rscript scripts/differential_expression.R
Rscript scripts/GO_enrichment.R
python scripts/correlation_analysis.py
```

Modify paths and parameters according to your local environment and dataset locations.

---

## Data Sources

The analyses are based on RNA-seq datasets from *S. solidus* across multiple developmental stages and hosts.

Relevant genomic resources include:

* [WormBase ParaSite](https://parasite.wormbase.org/Schistocephalus_solidus_prjeb527/Info/Index/?utm_source=chatgpt.com)
* [Transcriptome resource for Schistocephalus solidus](https://academic.oup.com/gigascience/article-abstract/doi/10.1186/s13742-016-0128-3/2720984?utm_source=chatgpt.com)

---


## Citation

If you use this repository or associated analyses, please cite the corresponding manuscript and repository.

Example:

```text
Gramolini L. et al.
Gene expression profiling across the three-host life cycle of Schistocephalus
solidus: how decoupled are the life stages?
```

---

## Contact

Maintainer: Laura Gramolini

GitHub: [LauraGramolini GitHub profile](https://github.com/LauraGramolini?utm_source=chatgpt.com)
