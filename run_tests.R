# 
# Author: lubo
# Created: Apr 18, 2018
###############################################################################

library(testthat)

devtools::load_all(".")

testthat::test_dir("tests/testthat")
