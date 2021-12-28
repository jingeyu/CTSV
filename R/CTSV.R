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
#' data(example_data)
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
#' 
#' # View on q-value matrix
#' head(result$qval)
#' 
#' # detect SV genes
#' re <- SVGene(result$qval,0.05)
#' #SV genes in each cell type:
#' print(re$SVGene)
#' 
#' @export
ctsv <- function(Y, loc, W, num_core=1){
    
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
        stop("Please remove cell types having zero proportion across all spots!")
    }
    if(nrow(Y)!=nrow(loc) | nrow(loc)!=nrow(W) | nrow(Y)!=nrow(W)){
        stop("Please keep the number of spots consistent in gene expression data matrix, location coordinate matrix and cell-type proportion matrix!")
    }
    if(sum(rownames(Y)!= rownames(loc))>0 | sum(rownames(W)!= rownames(loc))>0 | sum(rownames(Y)!= rownames(W))>0){
        stop("Please match spots' names in gene expression data matrix, location coordinate matrix, and cell-type proportion matrix!")
    }
    ## set number of core in run
    if(num_core > 1){
        if(num_core>detectCores()){warning("The number of cores you're setting is larger than detected cores!");
            num_core = detectCores()-1}
    }
    
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
    quan <- c(0.2,0.4,0.6,0.8,1.0)
    psi1 <- quantile(abs(S[,1]), quan)
    psi2 <- quantile(abs(S[,2]), quan)
    P_VAL <- array(NA, dim = c(G,2*K,5))
    
    cat(paste("## ======= CTSV input ====== \n"))
    cat(paste("## number of total spots: ", n,"\n"))
    cat(paste("## number of total genes: ", G,"\n"))
    
    cl<- makeCluster(num_core)
    registerDoSNOW(cl)
    pattern <- c("linear","gau1","gau2","cos1","cos2")
    for(fit_pat in pattern){
        if(fit_pat == "gau1"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
        }else if(fit_pat == "gau2"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[3]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[3]^2)
        }else if(fit_pat == "cos1"){
            h1 <- cos(2*pi*S[,1] / psi1[2])
            h2 <- cos(2*pi*S[,2] / psi2[2])
        } else if(fit_pat == "cos2"){
            h1 <- cos(2*pi*S[,1] / psi1[3])
            h2 <- cos(2*pi*S[,2] / psi2[3])
        }else{
            h1 <- S[,1]
            h2 <- S[,2]
        }
        cat(paste("## Fitting basis function: ", fit_pat,"\n"))
        
        Tmp <- cbind(W * h1, W * h2, W)
        colnames(Tmp) <- 1:ncol(Tmp)
        P_gene <- function(g,Y,Tmp,h1,h2){
            y <- Y[,g]
            K <- ncol(Tmp)/3
            ell <- rowSums(Y) / median(rowSums(Y))
            err <- try(fm_zinb0 <- zeroinfl(y ~ -1+offset(log(ell))+Tmp|1,
                                            dist = "negbin",link = "probit",
                                            control = zeroinfl.control(method = "CG"
                                                                       # reltol = 0.0001
                                                                       # maxit = 5000
                                            )), silent = TRUE)
            if(class(err) == 'try-error') {
                p_val <- rep(-1,2*K)
            } else{
                p_val <- coef(summary(fm_zinb0))$count[,4]
                nind <- 2*(length(p_val) - 1)/3
                p_val <- p_val[1:nind]
                
            }
            return(p_val)
            
        }
        
        P_VAL[,,match(fit_pat,pattern)] = foreach(g=1:G,
                                                  .combine=rbind,
                                                  .packages=c("pscl")) %dopar% P_gene(g,Y,Tmp,h1,h2)
        
        rownames(P_VAL[,,match(fit_pat,pattern)]) <- colnames(Y)
    }
    
    p_val <- foreach(g=1:G,
                     .combine=rbind,
                     .packages=c("pscl")) %dopar% P_gene(g,Y,Tmp,h1,h2)
    
    
    
    stopCluster(cl)
    # Cauchy combination rule
    P_VAL[which(is.na(P_VAL))] <- 1
    P_VAL[P_VAL == -1] <- 1
    P_VAL <- tan((0.5 - P_VAL)*pi)
    T_cau0 <- apply(P_VAL, c(1,2), mean)
    P_val <- 1-pcauchy(T_cau0)
    # convert q-values into q-values
    Q_val <- matrix(qvalue(c(P_val))$qvalue,G,2*K)
    rownames(Q_val) <- colnames(Y)
    
    return(list("pval" = P_val,
                "qval" = Q_val))
}


#======Get SV genes=========
#' @title Report spatially variable genes
#' @param Q_val: A G by 2K q-value matrix, where G is the number of genes and K is the number of cell types.
#' @param thre.alpha: numeric, a q-value threshold to control FDR less than thre.alpha.
#' @return A list with a G by 2K 0-1 matrix and a list with names of SV genes in each cell type. The first K columns of the 0-1 matrix correspond to the coordinate of \eqn{S_1}, and the last K columns to the coordinate of \eqn{S_2}.
#' \item{SV}{A G by 2K 0-1 matrix. The first K columns correspond to the coordinate of \eqn{S_1}, the last K columns to the coordinate of \eqn{S_2}.}
#' \item{SVGene}{A list with names of SV genes in each cell type.}
#' @export
SVGene <- function(Q_val, thre.alpha){
    G <- nrow(Q_val)
    K <- ncol(Q_val)/2
    SVmat <- matrix(0,G,2*K)
    SVmat[Q_val < thre.alpha] <- 1
    if(is.null(rownames(Q_val))){
        stop("Please name rows of q-value matrix with corresponding gene names!")
    }
    all_gene <- rownames(Q_val)
    svg <- list()
    for(k in 1:K){
       svg[[k]] <- all_gene[which(rowSums(SVmat[,c(k,k+K)]) > 0)]
    }
    return(list("SV" = SVmat,
                "SVGene" = svg))
}

