## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    crop = NULL 
)

## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE----------------
## Track time spent on making the vignette
startTime <- Sys.time()

## ----"install", eval = FALSE--------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#        install.packages("BiocManager")
#    }
#  
#  BiocManager::install("CTSV")
#  
#  ## Check that you have a valid Bioconductor installation
#  BiocManager::valid()

## ----library packages,results = 'hide'----------------------------------------
suppressPackageStartupMessages(library(CTSV))
suppressPackageStartupMessages(library(SpatialExperiment))

## ----ST-----------------------------------------------------------------------
data("CTSVexample_data", package="CTSV")
spe <- CTSVexample_data[[1]]
W <- CTSVexample_data[[2]]
gamma_true <- CTSVexample_data[[3]]
Y <- assay(spe)
# dimension of bulk ST data
dim(Y)
# dimension of cell-type proportion matrix:
dim(W)
# SV genes in each cell type:
colnames(Y)[which(gamma_true[,1] == 1)]
colnames(Y)[which(gamma_true[,2] == 1)]
# Number of SV genes at the aggregated level:
sum(rowSums(gamma_true)>0)

## ----Run CTSV-----------------------------------------------------------------
result <- CTSV(spe,W,num_core = 8)

## ----results------------------------------------------------------------------
# View on q-value matrix
head(result$qval)

## ----SVgenes------------------------------------------------------------------
re <- svGene(result$qval,0.05)
# SV genes in each cell type:
re$SVGene

## ----session information------------------------------------------------------
sessionInfo()

