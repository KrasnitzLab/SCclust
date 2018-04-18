# 
# Author: lubo
# Created: Apr 18, 2018
###############################################################################



library(futile.logger)
flog.threshold(DEBUG)

devtools::load_all(".")


options(testthat.output_file = "junit_test_results.xml")

to_run <- c(
    'tests/testthat/test_fileutils.R'
)

for(t in to_run) {
  testthat::test_file(t, reporter="junit")
}