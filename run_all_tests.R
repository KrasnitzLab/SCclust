library(devtools)
library(testthat)

devtools::load_all()

#testthat::test_file("tests/testthat/test_segment_nyu003.R")


#testthat::test_file("tests/testthat/test_fisher_fdr.R")
#testthat::test_file("tests/testthat/test_hclust_tree.R")
#testthat::test_file("tests/testthat/test_find_clones.R")
#testthat::test_file("tests/testthat/test_find_subclones.R")


#testthat::test_dir("tests/")

to_run <- c(
    'tests/testthat/test_segment_nyu003.R',
    'tests/testthat/test_chrom_numeric.R',
    'tests/testthat/test_find_clones.R',
    'tests/testthat/test_fisher_fdr.R',
    'tests/testthat/test_pinmat.R',
    'tests/testthat/test_utils.R',
    'tests/testthat/test_fileutils.R',
    'tests/testthat/test_find_subclones.R',
    'tests/testthat/test_hclust_tree.R',
    'tests/testthat/test_preprocess.R',
    'tests/testthat/test_simfisher.R')
#tests/testthat/test_simfisher_convergence_nyu003.R
#tests/testthat/test_simfisher_convergence_gl61.R

for(t in to_run) {
  testthat::test_file(t)
}
