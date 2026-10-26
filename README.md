# Gut Microbiome Biomarker Benchmark

A benchmark of computational methods for gut microbiome disease biomarker identification, accompanied by the exploration of best practices. 

This repository contains the R scripts and analysis pipelines for my undergraduate thesis: **"Comparison of computational methods for gut microbiome biomarker identification and exploration of best practices"** (Huazhong University of Science and Technology).

## Overview

Gut microbiome data are inherently sparse, compositional, and prone to batch confounders. These characteristics reduce the reproducibility of biomarker findings. This project constructs a unified and reproducible benchmarking pipeline using simulated data constrained by real-world features from the **Gut Microbiome Repository (GMrepo)**. 

We systematically evaluated 5 mainstream methods across 4 challenging scenarios:
*   **Methods:** ANCOM-BC2, DESeq2 (default/noCooks), LEfSe, MaAsLin2, and Random Forest (RF).
*   **Scenarios:** A (No confound), B (Batch confound), C (Extremely sparse/weak signal), D (Large effect).
*   **Metrics:** F1, FDR, Top-K hits, Jaccard stability, Batch error rate, and Runtime.

## Repository Structure

The scripts are organized by the analytical workflows of the thesis:

*   `mockdata_generate.R`: Simulates microbiome count data based on GMrepo statistical constraints (Dirichlet-Multinomial with zero-inflation and signal implantation).
*   **Method-specific Scripts:**
    *   `ancombc2.R`: Runs ANCOM-BC2 in both *naive* and *batch-adjusted* modes.
    *   `deseq2.R`: Runs DESeq2 (default and noCooks configurations).
    *   `lefse.R`: Runs LEfSe analysis.
    *   `maaslin2.R`: Runs MaAsLin2 with TSS normalization and log transformation.
    *   `rf.R`: Runs Random Forest with rigorous cross-validation and feature ranking.
*   `benchmark.R`: The core evaluation script computing metrics (Precision, Recall, F1, Jaccard, Top-K) and aggregating results across scenarios.
*   `FullScale_analysis`: Scripts for applying the recommended baseline pipeline to real GMrepo cohorts.
*   `representative_verification`: Scripts for the consistency testing and overlap analysis between different methods on real data.


## Data Availability
Due to the large file size and data usage policies of the GMrepo database, the raw count matrices and metadata are not included in this repository.
The scripts assume that the input data (simulated or real) are formatted as a taxa x sample count matrix and a corresponding metadata table. Please refer to the header comments in mockdata_generate.R for the expected data structure.


## Key Findings (Best Practices)
Based on our benchmark and real-cohort analysis, we propose the following workflow recommendations for GMrepo:

Default/Routine Analysis: Use MaAsLin2 (for scalability and covariate handling) and LEfSe (for non-parametric robustness) in conjunction. ANCOM-BC2 can serve as a supplementary method for compositional bias correction.

Sparse/Weak Signal Scenarios: Rely on Random Forest for Top-K candidate generation. Statistical methods may suffer from severe power loss.

Batch Confounding: If batch information is available, use covariate-adjusted models (e.g., MaAsLin2). If not, LEfSe exhibits lower batch mis-selection rates in naive analysis.

Output Strategy: Clearly distinguish between "Significance Sets" (for statistical evidence) and "Top-K Candidate Rankings" (for prioritization in downstream validation).


*Author: Zihan Sun (Undergraduate Thesis, HUST)*
