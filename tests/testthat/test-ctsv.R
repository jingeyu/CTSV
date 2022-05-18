context("correct input files for ctsv function")

data(CTSVexample_data)
Y <- CTSVexample_data[[4]]
W <- CTSVexample_data[[3]]
loc <- CTSVexample_data[[2]]

test_that("errors are returned when there are missing inputs",{
    
    expect_error(ctsv(),
                 "Please*")
    expect_error(ctsv(Y=Y),
                 "Please*")
    expect_error(ctsv(Y=Y,loc = loc),
                 "Please*")
    
})

test_that("errors are returned when expression data, loaction matrix, cell-type proportion matrix and number of cores are not in correct format",{
    
    expect_error(ctsv(Y = as.character(Y),loc = loc,W=W,num_core = 1),
                 "Please*")  
    expect_error(ctsv(Y = Y,loc = as.character(loc),W=W,num_core = 1),
                 "Please*")  
    expect_error(ctsv(Y = Y,loc = loc,W=as.character(W),num_core = 1),
                 "Please*")  
    expect_error(ctsv(Y = Y,loc = loc,W=W,num_core = 1.5),
                 "Please*")  
    expect_error(ctsv(Y = Y[-2,],loc = loc,W=W),
                 "Please*")  
    expect_error(ctsv(Y = Y,loc = loc[-3,],W=W),
                 "Please*")  
    expect_error(ctsv(Y = Y,loc = loc,W=W[-1,]),
                 "Please*")  
    spot <- rownames(Y)
    spot[1] <- "my"
    rownames(Y) <- spot
    expect_error(ctsv(Y = Y,loc = loc,W=W),
                 "Please*")  
})

