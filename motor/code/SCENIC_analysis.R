# Suppress loading messages when building the HTML
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
library(SCENIC)
org <- "mgi" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=1) 

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 