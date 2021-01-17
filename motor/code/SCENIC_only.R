# yum install xorg-x11-server-Xvfb
# xvfb-run R script

# SCENIC
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(KernSmooth)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
})

# Initialize SCENIC settings
org <- "mgi" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Read data & filter
dbs <- defaultDbNames[["mgi"]]
expMat = read.csv("data3/geneExp_motor.csv")
NMA_result = read.csv("data3/efeatures_NMA.csv")
expMat = expMat[!duplicated(expMat[,1]),]
rownames(expMat) = expMat[,1];expMat =as.matrix(expMat[,-1])
expMat = expMat[,NMA_result$cellnames[NMA_result$ttype == "Pvalb"]]

# Gene/cell filter/selection
ne_genes = read.csv("data3/neuronal_genes.csv")
TF_genes = read.csv("data3/TFgenes_simp.csv")
genesKept = union(ne_genes[,1],TF_genes[,1])
genesKept = intersect(genesKept,rownames(expMat))
expMat_filtered <- expMat[genesKept, ]
dim(expMat_filtered)
expMat_filtered <- log2(expMat_filtered+1)

# Correlation
runCorrelation(expMat_filtered, scenicOptions)
# GENIE3
runGenie3(expMat_filtered, scenicOptions)

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

options(bitmapType = 'cairo')
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
runSCENIC_3_scoreCells(scenicOptions, expMat_filtered)
#saveRDS(scenicOptions, file="int/scenicOptions.Rds")