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

- `igraph` (R package) https://igraph.org/r/ or `cytoscape` (open source platform) https://cytoscape.org/ to draw igraphs for cluster specific gene regulatory networks.

- `g:Profiler` A web server for functional enrichment analysis. Available at https://biit.cs.ut.ee/gprofiler/gost.

## Demo for aligning single cell multi-modal data

A demo for aligning single-cell multi-modal data is available at `demo_112_visual_cortex/` for a Patch-seq dataset (112 cells) in the mouse visual cortex https://github.com/berenslab/layer4. The input data includes the gene expression levels and electrophysiological features of these cells in a rda file. Also, we provide a list of neuronal marker genes as select gene features. The code performs dimensional reduction on single-cell electrophysiological data and gene expression with three major methods: Nonlinear manifold alignment (NMA), Canonical correspondence analysis (CCA) and Principal component analysis (PCA).  The expected output includes the visualization of the latent space after each alignment method. The total running time of this demo was less than 5 minutes on a local laptop.

## Instructions for use for other datasets

### 3654 cells in the mouse visual cortex

The data and codes for analyzing a large single-cell multi-modal data in the mouse visual cortex analysis (N=3654 cells) is at `mouse_visual_cortex/`. The raw gene expression data, `20200513_Mouse_PatchSeq_Release_count.v2.csv.tar`, along with metadata is available at http://data.nemoarchive.org/other/grant/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/ . Also, the raw electrophysiological data can be obtained through the following command.

```python
pip install dandi
dandi download https://dandiarchive.org/dandiset/000020
```

After downloading the raw electrophysiological data (around 112 GB), we can extract the electrophysiological features by the python code, `extracting_e_features_visual.py`. Alternatively, you can directly use the processed data, `visual_data_filtered.rda by data_clean_visual.R`, which include processed gene expression, select electrophysiological features (e-features) and cell type information for all 3654 cells.

The codes for nonlinear manifold alignment (NMA) and clustering is available in `dimreduction_clustering_visual.Rmd`, which outputs two tables representing the NMA's latent spaces of genes and e-features. The codes for identifying the differential expressed genes and cluster-specific e-features are also available in `dimreduction_clustering_visual.Rmd`. The code file, `regulatory_network_visual.R` predicts gene regulatory network for cross-modal cell clusters, and should be run on linux. The code file, `predicting_efeatures_visual.rmd` predicts the cluster specific e-features by the method for predicting differential expressed genes, and also generates a plot showing the $R^2$ for e-features prediction models.

Others plots can be generated through the following codes files: `igraph_visual.Rmd` visualizes the gene regulatory network by igraph for each cluster; `bibiplot_visual.ipynb` and `regression_on_latent_space_CV_visual.ipynb` generate the bibiplots and fit the regression models for the gene expression and e-feature latent spaces; `network_plot_visual.Rmd` visualizes the expression relationships and calculates the expression correlations of three marker genes in cluster 1: "Irf5","Spi1" and "Atf3"; `mtype_plot_visual.r` plots the m-type trajectory on the NMA's latent space after adding 142 spiny (excitatory) cells into the NMA alignment.

### 1208 cells in the mouse motor cortex

The single-cell multi-modal data in the mouse motor cortex (N=1208 cells) was obtained from https://github.com/berenslab/mini-atlas. We also include them in `mouse_motor_cortex/data/`. The codes for generating our results for this dataset are available in `mouse_motor_cortex/codes`, and the code file names follow the same naming strategy with the mouse visual cortex but replace the suffix by '_motor'.

