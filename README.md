# RBP Network Analysis Framework

Repository containing R scripts for validation, benchmarking, prioritization, enrichment analysis, and visualization of RNA-binding protein (RBP) regulatory interaction networks derived from ENCODE perturbation datasets.

---

## Overview

This repository provides a downstream analysis framework for studying RBP regulatory interaction networks using multiple complementary validation and benchmarking approaches.

The workflow integrates:

- consensus network construction
- co-expression validation
- eCLIP-supported interaction analysis
- disease enrichment analysis
- topology benchmarking
- cross-cell-line comparison
- synthetic dataset benchmarking
- interaction prioritization
- publication-ready visualization

The scripts were cleaned and standardized for publication and reproducible use.

---

# Main Scripts

## 01_data_preprocessing_and_matrix_generation.R
Preprocessing and preparation of ENCODE-derived perturbation transcriptomic datasets.

## 02_network_inference_methods.R
Application of multiple computational network inference approaches.

## 03_consensus_network_construction.R
Construction of consensus regulatory networks from multiple inference methods.

## 04_external_validation_and_crosscomparison.R
Cross-validation of inferred interactions using external resources.

## 05_coexpression_validation_plots.R
Co-expression validation using TCGA LIHC and GTEx liver datasets.

## 06_network_module_enrichment_analysis.R
Community detection and functional enrichment analysis of RBP interaction modules.

## 07_network_topology_comparison.R
Benchmarking of inferred network topology against reference regulatory networks.

## 08_eclip_target_network_construction.R
Construction of eCLIP-supported regulator-target interaction networks.

## 09_link_prioritization_validation_scoring.R
Multi-layer validation and prioritization of candidate regulatory interactions.

## 10_supplementary_figure_network_benchmarking.R
Generation of supplementary benchmarking figures using simulated ENCODE-like datasets.

## 11_ranked_interaction_disease_enrichment_analysis.R
Disease enrichment analysis of ranked interactions using threshold and sliding-window approaches.

## 12_archived_cross_cellline_consensus_network_analysis.R
Archived exploratory workflow for cross-cell-line consensus network comparison.

## 13_syncode_synthetic_dataset_similarity_analysis.R
Benchmarking of SYNCODE-generated synthetic transcriptomic datasets against real perturbation datasets.

## 14_multifeature_interaction_scoring_and_visualization.R
Integration, scoring, ranking, visualization, and enrichment analysis of regulatory interactions using multiple evidence layers.

---

# Utility Functions

The repository also contains a collection of helper functions used throughout the workflow for:

- overlap statistics
- consensus network construction
- network comparison
- auxiliary downstream analyses

---

# Main Dependencies

Core R packages used throughout the project include:

- tidyverse
- data.table
- igraph
- ggplot2
- ggpubr
- caret
- rpart
- mclust
- enrichplot
- DOSE
- org.Hs.eg.db
- RCAS
- GenomicRanges
- rtracklayer
- pheatmap
- R.matlab
- gprofiler2

---

# Input Data

The workflow relies on multiple external datasets, including:

- ENCODE perturbation transcriptomics
- ENCODE eCLIP datasets
- GTEx liver expression data
- TCGA LIHC expression and survival data
- GENCODE annotations
- TRRUST regulatory networks
- FunCoup interaction networks
- MSigDB functional gene sets (selected)

---

# Expected Outputs

The workflow generates:

- inferred regulatory networks
- consensus adjacency matrices
- ranked interaction tables
- disease enrichment summaries
- benchmarking figures
- Cytoscape-compatible network exports
- publication-ready figures
- supplementary figures

---

# Citation

If you use this repository in your work, please cite the associated publication and/or Zenodo archive.
