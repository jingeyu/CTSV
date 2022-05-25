\name{CTSVexample_data}
\alias{CTSVexample_data}
\docType{data}
\title{
A simulated data set
}   

\format{A list containing 3 elements.}

\description{A simulated data set for demonstrating how to use the ctsv function.
\describe{
  \item{gamma_true}{An 20 by 2 0-1 matrix, indicator of SV genes.}
  \item{spe}{A SpatialExperiment class.}
  \item{W}{An 100 by 2 cell-type proportion matrix.}
}
}

\examples{
\donttest{
library(SpatialExperiment)
rDirichlet <- function(n,alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
    
}
seed <- 20210509
set.seed(seed)
#  gene numbers
G <- 20
# cell type numbers
K <- 2
# spot numbers
n <- 100
# number of DE genes
DE_num <- 10
# drop out probability
pai <- 0.5
# parameter of NB distribution
size = 100

# coordinates of spots
loc <- NULL
for(i in 1:10){
    for(j in 1:10){
        loc <- rbind(loc,c(i,j))
    }
}
rownames(loc) <- paste0("spot",1:n)
colnames(loc) <- c("x","y")
NDE_scrna <- rnorm(G, mean=2, sd=0.2)
scrna_1 <- NDE_scrna
scrna_2 <- NDE_scrna
scrna_2[sample(1:G,DE_num,replace = FALSE)] <- rnorm(DE_num, mean=3, sd=0.2)
eta <- cbind(scrna_1,scrna_2)

gamma_true <- matrix(0, G, K)
gamma_true[11:13,1] <- 1
gamma_true[14:16,2] <- 1
beta1 <- matrix(0, G, K)
beta2 <- matrix(0, G, K)

# cell type proportion
W <- rDirichlet(n, c(1,2))
W <- t(W)

S <- t(loc) - colMeans(loc)
S <- t(S / apply(S, 1, sd))

h1 <- S[,1]
h2 <- S[,2]
beta1[gamma_true == 1] <- 1
beta2[gamma_true == 1] <- 0.5

log_lambda <- eta %*% W + beta1 %*% t(t(W) * h1) + beta2 %*% t(t(W) * h2)
W <- t(W)
Y <- matrix(rnbinom(G*n,size = size, mu = exp(c(log_lambda))), G, n)
set.seed(5)
r_unif <- matrix(runif(G*n),G,n)
Y[r_unif <= pai] <- 0
colnames(Y) = rownames(loc)
rownames(W) = rownames(loc)
rownames(Y) <- paste0("gene",1:G)
spe <- SpatialExperiment(
    assay = list(counts = Y),
    colData = loc,
    spatialCoordsNames = c("x","y")
)
CTSVexample_data <- list(spe,W,gamma_true)
}
}