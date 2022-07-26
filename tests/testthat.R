library(devtools)
library(testthat)
library(futile.logger)

flog.threshold(INFO)

devtools::load_all(".")
test_dir("tests/testthat", reporter=SummaryReporter$new())
