################################################################
# Internal functions to be used in the CTSV function
################################################################
.P_gene <- function(g,Y,Tmp,h1,h2){
    y <- Y[,g]
    K <- ncol(Tmp)/3
    ell <- rowSums(Y) / median(rowSums(Y))
    err <- try(fm_zinb0 <- pscl::zeroinfl(y ~ -1+offset(log(ell))+Tmp|1,
                                          dist = "negbin",link = "probit",
                                          control = zeroinfl.control(method = "CG"
                                          )), silent = TRUE)
    is_err <- methods::is(err,"class")
    if(is_err){
        p_val <- rep(-1,2*K)
    }else{
        p_val <- coef(summary(fm_zinb0))$count[,4]
        p_val <- p_val[seq_len(2*K)]
    }
    return(p_val)
    
}

#' @title Detection of cell-type-specific spatially variable genes
#' @param Y An n by G bulk ST raw count data matrix, where n is the spot number and G is the gene number. Each row stands for a spot and each column represents a gene. The column names of Y (i.e. gene names) must be specified, and the row names of Y must be consistent with row names of loc and W.
#' @param loc An n by 2 location matrix, where each row corresponds to a spot and records the two-dimensional coordinate of that spot.rd
#' @param W An n by K cell-type proportion matrix, where K is the number of cell types. The column names of W are cell type names.
#' @param num_core Number of cores if using paralleling. The default is one.
#' @return A list with a G by 2K matrix of p-values and a G by 2K matrix of q-values.
#' \item{pval}{A G by 2K matrix of p-values. The first K columns correspond to the first coordinate, and the last K columns to the second coordinate.}
#' \item{qval}{A G by 2K matrix of q-values. The first K columns correspond to the first coordinate, and the last K columns to the second coordinate.}
#' @examples
#' library(CTSV)
#' #read example data
#' data(CTSVexample_data)
#' Y <- CTSVexample_data[[4]]
#' W <- CTSVexample_data[[3]]
#' loc <- CTSVexample_data[[2]]
#' gamma_true <- CTSVexample_data[[1]]
#' # gene number
#' G <- ncol(Y)
#' # spot number
#' n <- nrow(Y)
#' # cell type number
#' K <- ncol(W)
#' 
#' # SV genes in each cell type:
#' print(colnames(Y)[which(gamma_true[,1] == 1)])
#' print(colnames(Y)[which(gamma_true[,2] == 1)])
#' # Number of SV genes at the aggregated level:
#' print(sum(rowSums(gamma_true)>0))
#' 
#' #--- Run CTSV ----
#' result <- ctsv(Y,loc,W,num_core = 8)
#' # View on q-value matrix
#' head(result$qval)
#' # detect SV genes
#' re <- svGene(result$qval,0.05)
#' #SV genes in each cell type:
#' print(re$SVGene)
#' @export
ctsv <- function(Y, loc, W, num_core=1){
    if (missing(Y)) {
        stop("Please include gene expression data!")
    }
    if (missing(loc)) {
        stop("Please include location matrix!")
    }
    if (missing(W)) {
        stop("Please include cell-type proportion matrix!")
    }
    if(as.integer(num_core)!=as.numeric(num_core)){
        stop("Please input integer num of cores!")
    }
    if(!is.matrix(Y) | !is.matrix(W) | !is.matrix(loc)){
        stop("Please input the matrix type of \"Y\", \"W\", and \"loc\"!")
    }
    if(is.null(colnames(Y))){
        stop("Please give colnames for the spatial expression matrix (gene names)!")
    }
    if(sum(is.na(Y))>0 | sum(is.na(loc))>0 | sum(is.na(W)) > 0){
        stop("Please remove NaNs in the datasets!")
    }
    if(sum(rowSums(Y) == 0)>0){
        stop("Please remove genes having zero values across all spots!")
    }
    if(sum(colSums(Y) == 0)>0){
        stop("Please remove spots having zero values across all genes!")
    }
    if(sum(colSums(W) == 0)>0){
        stop("Please remove cell types having zero proportion across all spots!")
    }
    if(sum(rowSums(W) == 0)>0){
        stop("Please remove spots having zero proportion across all cell types!")
    }
    if(nrow(Y)!=nrow(loc) | nrow(loc)!=nrow(W) | nrow(Y)!=nrow(W)){
        stop("Please keep the number of spots consistent in gene expression data matrix, location coordinate matrix and cell-type proportion matrix!")
    }
    if(sum(rownames(Y)!= rownames(loc))>0 | sum(rownames(W)!= rownames(loc))>0 | sum(rownames(Y)!= rownames(W))>0){
        stop("Please match spots' names in gene expression data matrix, location coordinate matrix, and cell-type proportion matrix!")
    }
    # make sure the sum of cell type proportions is equal to 1 in each spot.    
    W <- W / rowSums(W)
    # number of genes
    G <- ncol(Y)
    # number of spots
    n <- nrow(loc)
    # number of cell types
    K <- ncol(W)
    # normalize cell-type proportion matrix W to ensure the summation across cell types in one spot is equal to one.
    W <- W / rowSums(W)
    # Center and normalize coordinates of spots to have mean zero and standard deviation one.
    S <- t(loc) - colMeans(loc)
    S <- t(S / apply(S, 1, sd))
    quan <- c(0.4,0.6)
    psi1 <- quantile(abs(S[,1]), quan)
    psi2 <- quantile(abs(S[,2]), quan)
    P_VAL <- array(NA, dim = c(G, 2*K, 5))
    bp_param <- BiocParallel::MulticoreParam(workers=num_core)
    pattern <- c("linear","gau1","gau2","cos1","cos2")
    for(fit_pat in pattern){
        if(fit_pat == "gau1"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[1]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[1]^2)
        }else if(fit_pat == "gau2"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
        }else if(fit_pat == "cos1"){
            h1 <- cos(2*pi*S[,1] / psi1[1])
            h2 <- cos(2*pi*S[,2] / psi2[1])
        } else if(fit_pat == "cos2"){
            h1 <- cos(2*pi*S[,1] / psi1[2])
            h2 <- cos(2*pi*S[,2] / psi2[2])
        }else{
            h1 <- S[,1]
            h2 <- S[,2]
        }
        Tmp <- cbind(W * h1, W * h2, W)
        colnames(Tmp) <- seq_len(ncol(Tmp))
        res <- do.call(rbind,BiocParallel::bplapply(seq_len(G),.P_gene,BPPARAM = bp_param,Y=Y,Tmp = Tmp,h1=h1,h2=h2))
        P_VAL[,,match(fit_pat,pattern)] <- res
        rownames(P_VAL[,,match(fit_pat,pattern)]) <- colnames(Y)
    }
    P_VAL[which(is.na(P_VAL))] <- 1
    P_VAL[P_VAL == -1] <- 1
    P_VAL <- tan((0.5 - P_VAL)*pi)
    T_cau0 <- apply(P_VAL, c(1,2), mean)
    P_val <- 1-pcauchy(T_cau0)
    # convert q-values into q-values
    Q_val <- matrix(qvalue(c(P_val))$qvalue, G, 2*K)
    rownames(Q_val) <- colnames(Y)
    return(list("pval" = P_val,
                "qval" = Q_val))
}


#======Get SV genes=========
#' @title Report spatially variable genes
#' @param Q_val A G by 2K q-value matrix, where G is the number of genes and K is the number of cell types.
#' @param thre.alpha numeric, a q-value threshold to control FDR less than thre.alpha.
#' @return A list with a G by 2K 0-1 matrix and a list with names of SV genes in each cell type. The first K columns of the 0-1 matrix correspond to the coordinate of \eqn{S_1}, and the last K columns to the coordinate of \eqn{S_2}.
#' \item{SV}{A G by 2K 0-1 matrix. The first K columns correspond to the coordinate of \eqn{S_1}, the last K columns to the coordinate of \eqn{S_2}.}
#' \item{SVGene}{A list with names of SV genes in each cell type.}
#' @examples 
#' library(CTSV)
#' # Simulate a Q value matrix
#' K <- 2 # cell-type number
#' G <- 10 # gene number
#' set.seed(1)
#' Q_val <-matrix(runif(G*K,0,0.1),G,K)
#' rownames(Q_val) <- paste0("gene",seq_len(G))
#' # detect SV genes
#' re <- svGene(Q_val,0.05)
#' #SV genes in each cell type:
#' print(re$SVGene)
#' @export
svGene <- function(Q_val, thre.alpha=0.05){
    if (missing(Q_val)) {
        stop("Please include gene expression data!")
    }
    if(is.null(rownames(Q_val))){
        stop("Please name rows of q-value matrix with corresponding gene names!")
    }
    if(thre.alpha < 0 | thre.alpha > 1){
        stop("Please limit the threshold between the rangeo of 0 and 1.")
    }
    G <- nrow(Q_val)
    K <- ncol(Q_val)/2
    SVmat <- matrix(0, G, 2*K)
    SVmat[Q_val < thre.alpha] <- 1
    
    all_gene <- rownames(Q_val)
    svg <- list()
    for(k in seq_len(K)){
        svg[[k]] <- all_gene[which(rowSums(SVmat[, c(k, k+K)]) > 0)]
    }
    return(list("SV" = SVmat,
                "SVGene" = svg))
}

