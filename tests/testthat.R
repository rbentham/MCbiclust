Sys.setenv("R_TESTS" = "")

library(testthat)
library(MCbiclust)

test_check("MCbiclust")
