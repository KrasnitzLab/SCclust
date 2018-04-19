# 
# Author: lubo
# Created: Apr 18, 2018
###############################################################################

library(testthat)

library(futile.logger)
flog.threshold(DEBUG)

devtools::load_all(".")


options(testthat.output_file = "junit_test_results.xml")
testthat::test_dir("tests/testthat", reporter="junit")
