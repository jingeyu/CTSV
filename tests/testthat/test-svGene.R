context("correct input files for svGene function")

set.seed(1)
Q_val <-matrix(runif(20,0,0.1),10,2)
rownames(Q_val) <- paste0("gene",seq_len(10))
test_that("errors are returned when there are missing inputs",{
    expect_error(svGene(),
                 "Include*")
})

test_that("errors are returned when the q-value matrix and threshold are not in correct format",{
    expect_error(svGene(Q_val=Q_val, thre.alpha=2),
                 "The threshold*")  
    result <- svGene(Q_val=Q_val,0.05)
    expect_equal(is.list(result), TRUE)
    expect_equal(length(result),2)
    expect_equal(dim(Q_val),dim(result$SV))
    rownames(Q_val) <- NULL
    expect_error(svGene(Q_val=Q_val),
                 "Name*")  
})
