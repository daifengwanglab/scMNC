# Manifold learning analysis suggests novel strategies for aligning single-cell multi-modalities and revealing functional genomics for cellular electrophysiology

## Summary
Recent single-cell multi-modal data reveal different characteristics of single cells, such as transcriptomics, morphology, and electrophysiology. However, our understanding of functional genomics and gene regulation leading to various cellular characteristics remains elusive. To address this, we applied multiple machine learning approaches to align gene expression and electrophysiological data of single neuronal cells in the mouse visual and motor cortices. We found that nonlinear manifold learning outperforms others. After manifold alignment, the cell clusters highly correspond to transcriptomic and morphological cell-types, suggesting a strong nonlinear linkage between gene expression and electrophysiology at the cell-type level. Also, the aligned cells form developmental trajectories, showing the underlying neuronal development and continuous characteristic changes. We further clustered aligned cells into multiple cross-modal cell clusters and found that the clusters' differentially expressed genes can predict many electrophysiological features. Additional functional enrichment and gene regulatory network analyses for our cell clusters revealed potential novel molecular mechanistic insights from genes to electrophysiology at cellular resolution. 

## Flow chart
![alt text](https://github.com/daifengwanglab/scMNC/blob/main/cover_figure.png)

## System requirements

The analysis is based on R 4.02 and python 3.7. For the Single-Cell manifold alignment, GMM clustering and most follow up analysis, you will only need a standard computer with enough RAM to support the operations. For the gene regulatory Network, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Installation guide

Several packages are needed for the whole analysis.

- `ipfx 1.0.2` A python package designed for Allen SDK Neurodata Without Borders 2.0 (NWB) files to extract electrophysiology data files for the PatchSeq experiments, data access and installation gudiance is avaiable at [https://ipfx.readthedocs.io/en/latest/tutorial.html].


