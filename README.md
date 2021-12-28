## CTSV: Identification of cell-type-specific spatially variable genes accounting for excess zeros

The R package CTSV implements the CTSV approach developed by Jinge Yu and Xiangyu Luo that detects cell-type-specific spatially variable genes accounting for excess zeros. CTSV directly models sparse raw count data through a zero-inflated negative binomial regression model, incorporates cell-type proportions, and performs hypothesis testing based on R package pscl. The package outputs p-values and q-values for genes in each cell type, and CTSV is scalable to datasets with tens of thousands of genes measured on hundreds of spots. CTSV can be installed in Windows, Linux, and Mac OS. 


## Prerequisites and Installation

1. R version >= 3.6.
2. R packages:

  CRAN package: stats (>= 4.1.0), doParallel (>= 1.0.16), doSNOW (>= 1.0.19), foreach (>= 1.5.1), pscl (>= 1.5.5)

  Bioconductor package: qvalue (>=2.24.0)
  
3. Install the package CTSV.

```
devtools::install_github("jingeyu/CTSV")
```


## Example Code
The following shows an example that runs the main functions "ctsv" and "SVGene" in our package. 

``` {r, eval=FALSE}
library(CTSV)

#read example data
data(example_data)
# gene number
G <- ncol(Y)
# spot number
n <- nrow(Y)
# cell type number
K <- ncol(W)

# SV genes in each cell type:
print(colnames(Y)[which(gamma_true[,1] == 1)])
print(colnames(Y)[which(gamma_true[,2] == 1)])
# Number of SV genes at the aggregated level:
print(sum(rowSums(gamma_true)>0))

#--- Run CTSV ----
result <- ctsv(Y,loc,W,num_core = 8)

# View on q-value matrix
head(result$qval)

# detect SV genes
re <- SVGene(result$qval,0.05)
#SV genes in each cell type:
print(re$SVGene)
```
or you can simply run
``` {r, eval=FALSE}
library(CTSV)
example(ctsv)
```

## Remarks
* If you have any questions regarding this package, please contact Jinge Yu at yjgruc@ruc.edu.cn.

