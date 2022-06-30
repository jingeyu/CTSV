context("correct input files for ctsv function")
data(CTSVexample_data)
spe <- CTSVexample_data[[1]]
W <- CTSVexample_data[[2]]
gamma_true <- CTSVexample_data[[3]]
test_that("errors are returned when there are missing inputs",{
    expect_error(CTSV(),
                 "Include*")
    expect_error(CTSV(spe),
                 "Include*")
})  


test_that("errors are returned when expression data, loaction matrix, cell-type proportion matrix and number of cores are not in correct format",{
    
    expect_error(CTSV(spe = W,W=W,num_core = 1),
                 "Include*")  
    expect_error(CTSV(spe = spe,W=as.character(W),num_core = 1),
                 "Include*")  
    expect_error(CTSV(spe = spe,W=W,num_core = 1.5),
                 "Input*")  
    expect_equal(ncol(spe), nrow(W))

})

test_that("The return of CTSV function",{
    result <- CTSV::CTSV(spe = spe, W = W, num_core = 8)
    expect_equal(is.list(result),
                 TRUE)
    expect_equal(length(result),2)
    expect_equal(ncol(result$pval),2*ncol(W))
    expect_equal(nrow(result$pval),nrow(spe))
    expect_equal(dim(result$pval), dim(result$qval))
})  



test_that("The input problems of CTSV function",{
spot <- rownames(W)
spot[1] <- "my"
rownames(W) <- spot
expect_error(CTSV(spe = spe,W=W),
             "Keep*")  
})