## ----setup,results = 'hide'---------------------------------------------------
library(CTSV,quietly = TRUE)
library(doSNOW)
library(doParallel)
library(foreach)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  out.width = "100%"
)

## ----ST-----------------------------------------------------------------------
dir <- system.file(package = 'CTSV') #directory for the example data
load(paste0(dir,"/example_data.RData"))
# zero rate in data:
print(summary(rowMeans(Y==0)))
# dimension of bulk ST data
print(dim(Y))
# dimension of cell-type proportion matrix:
print(dim(W))
# SV genes in each cell type:
print(colnames(Y)[which(gamma_true[,1] == 1)])
print(colnames(Y)[which(gamma_true[,2] == 1)])
# Number of SV genes at the aggregated level:
print(sum(rowSums(gamma_true)>0))

## ----Run CTSV-----------------------------------------------------------------
result <- ctsv(Y,loc,W,num_core = 1)

## ----results------------------------------------------------------------------
# View on q-value matrix
head(result$qval)

## ----SVgenes------------------------------------------------------------------
re <- SVGene(result$qval,0.05)
# SV genes in each cell type:
print(re$SVGene)

