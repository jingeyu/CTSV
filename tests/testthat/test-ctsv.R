context("correct input files for ctsv function")
data(CTSVexample_data)
spe <- CTSVexample_data[[1]]
W <- CTSVexample_data[[2]]
gamma_true <- CTSVexample_data[[3]]
test_that("errors are returned when there are missing inputs",{
    expect_error(ctsv(),
                 "Please*")
    expect_error(ctsv(spe),
                 "Please*")
})  

test_that("errors are returned when expression data, loaction matrix, cell-type proportion matrix and number of cores are not in correct format",{
    
    expect_error(ctsv(spe = W,W=W,num_core = 1),
                 "Please*")  
    expect_error(ctsv(spe = spe,W=as.character(W),num_core = 1),
                 "Please*")  
    expect_error(ctsv(spe = spe,W=W,num_core = 1.5),
                 "Please*")  
    spot <- rownames(W)
    spot[1] <- "my"
    rownames(W) <- spot
    expect_error(ctsv(spe = spe,W=W),
                 "Please*")  
})

