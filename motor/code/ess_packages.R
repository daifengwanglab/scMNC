if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
BiocManager::install("SingleCellExperiment")
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
#BiocManager::install(c("AUCell", "RcisTarget"))
#BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/AUCell")
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/GENIE3")
devtools::install_github("aertslab/SCENIC@v1.1.2")
packageVersion("SCENIC")

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))


# To export/visualize in http://scope.aertslab.org
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
