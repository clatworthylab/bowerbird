library(testthat)
library(bowerbird)
library(reticulate)
reticulate::use_condaenv("r-reticulate", required = TRUE)

test_check("bowerbird")
