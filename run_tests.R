# 
# Author: lubo
# Created: Apr 18, 2018
###############################################################################


library(devtools)
library(testthat)
library(futile.logger)

devtools::load_all()

flog.threshold(DEBUG)

options(testthat.output_file = "junit_test_results.xml")

testthat::test_dir("tests/", reporter="junit")

