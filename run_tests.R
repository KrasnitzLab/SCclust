library(devtools)
library(testthat)

devtools::load_all()

# testthat::test_file("tests/testthat/test_fisher_fdr.R")
testthat::test_file("tests/testthat/test_hclust_tree.R")

# testthat::test_dir("tests/")
