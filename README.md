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

## Mouse Visual Cortex

The data and code for visual cortex analysis is at `mouse_visual_cortex/`. The original data is too large in size that can not be uploaded to github. To get the original data, you can go to http://data.nemoarchive.org/other/grant/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/ to get the gene expression and metadata, among which the gene expression should be `20200513_Mouse_PatchSeq_Release_count.v2.csv.tar`. You can get the electrophysiological data through the following command. 

```python
pip install dandi
dandi download https://dandiarchive.org/dandiset/000020
```

After downloading the raw data (around 112 GB), you can extract the electrophysiological features by `extracting_e_features_visual.py`. Alternatively, you can just use the cleaned data by `data_clean_visual.R` that stored in `visual_data_filtered.rda`, which contains gene expression, electrophysiological features and cell type labels for 3654 cells.

The nonlinear manifold alignment and clustering procedure is coded in `dimreduction_clustering_visual.Rmd`, which outputs two tables representating two features' latent space. The code also renders the differential expressed genes and cluster specific e-features. `regulatory_network_visual.R` Get the TF genes and gene regulatory network for specific clusters, note that this code should be run on a *linux* system. `predicting_efeatures_visual.rmd` predict the cluster specific e-features through differential expressed genes, it will render a plot containing the $R^2$ for e-features prediction models.

Several plots are made through the following code: `igraph_visual.Rmd` gets the igraph for network in each cluster; `bibiplot_visual.ipynb` and `regression_on_latent_space_CV_visual.ipynb` draw the bibiplot and fit regression models on gene expression and e-feature latent space; `network_plot_visual.Rmd` draws the correlation of a special case of three marker genes in cluster 1: "Irf5","Spi1" and "Atf3"; `mtype_plot_visual.r` draws the m-type trajectory on latent space when adding 142 spiny (excitatory) cells into analysis.

## Mouse Visual Cortex

The motor cortex data is obtained from https://github.com/berenslab/mini-atlas. `m1_patchseq_exon_counts.csv` and `m1_patchseq_intron_counts.csv` contains the gene expression data, `m1_patchseq_ephys_features.csv` contains the e-features, and `m1_patchseq_meta_data.csv` is the meta data for cell types. The naming rule for motor code is similar to the visual cortex, so it would be easy to find what problem each code is dealing with.




