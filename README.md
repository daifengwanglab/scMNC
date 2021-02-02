# Manifold learning analysis suggests novel strategies for aligning single-cell multi-modalities and revealing functional genomics for cellular electrophysiology

## Summary
Recent single-cell multi-modal data reveal different characteristics of single cells, such as transcriptomics, morphology, and electrophysiology. However, our understanding of functional genomics and gene regulation leading to various cellular characteristics remains elusive. To address this, we applied multiple machine learning approaches to align gene expression and electrophysiological data of single neuronal cells in the mouse visual and motor cortices. We found that nonlinear manifold learning outperforms others. After manifold alignment, the cell clusters highly correspond to transcriptomic and morphological cell-types, suggesting a strong nonlinear linkage between gene expression and electrophysiology at the cell-type level. Also, the aligned cells form developmental trajectories, showing the underlying neuronal development and continuous characteristic changes. We further clustered aligned cells into multiple cross-modal cell clusters and found that the clusters' differentially expressed genes can predict many electrophysiological features. Additional functional enrichment and gene regulatory network analyses for our cell clusters revealed potential novel molecular mechanistic insights from genes to electrophysiology at cellular resolution. 

## Flow chart
![alt text](https://github.com/daifengwanglab/scMNC/blob/main/cover_figure.png)

## System requirements

The analysis is based on R 4.02, python 3.7 and several web based softwares. For the Single-Cell manifold alignment, GMM clustering and most follow up analysis, you will only need a standard computer with enough RAM to support the operations. For the gene regulatory network, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Package and software guide

Several packages or softwares are needed for the whole analysis.

- `ipfx 1.0.2` A python package designed for Allen SDK Neurodata Without Borders 2.0 (NWB) files to extract electrophysiological data files for the PatchSeq experiments, data access and installation gudiance is available at https://ipfx.readthedocs.io/en/latest/tutorial.html.

- `ManiNetCluster` An R package depending on depends on numpy, scipy, matplotlib, and scikit-learn. To install, a computer will need both the R statistical computing environment (version 3.4.1 or later) and Python (version 3.6.1 or later). We will need the mainfunction ManinetCluster in this package to perform multiple manifold learning dimensional reduction methods. Installation gudiance is available at https://github.com/daifengwanglab/ManiNetCluster.

- `Seurat` An R package designed for QC, analysis, and exploration of single-cell RNA-seq data. We will utlize the `FindAllMarkers` to find gene markers and electrophysiological markers, introduction could be found at https://satijalab.org/seurat/. 

- `SCENIC` An R package for Single-Cell Regulatory Network Inference, installation gudiance is available at https://github.com/aertslab/SCENIC.

- `igraph` (R package) or `cytoscape` (open source platform) to draw igraphs for cluster specific gene regulatory networks.

- `g:Profiler` A web server for functional enrichment analysis. Available at https://biit.cs.ut.ee/gprofiler/gost.

## 112 Demo data: to begin

You can get familiar with our analysis process through the 112 single cell patch-seq data from visual cortex, which could be accessed through `demo_102_visual_cortex/code/`. This code performs dimensional reduction on single-cell electrophysiological data and gene expression with 3 methods: Nonlinear manifold alignment (NMA), Canonical correspondence analysis (CCA) and Principal component analysis (PCA). The total running time of this demo should be less than 5 minutes.




