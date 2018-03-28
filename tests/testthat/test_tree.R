# 
# Author: lubo
# Created: Mar 28, 2018
###############################################################################


context("check tree_py output")

test_that("check tree_py works as expected", {
  data_dir <- Sys.getenv("SGAINS_DATA")

  mdist_filename <- file.path(
      data_dir,
      "gleason6.1/results/",
      "GL6.1smear1bpLog10FisherP.txt")
  expect_true(file.exists(mdist_filename))
  mdist <- load_table(mdist_filename)
  
  
  # res <- TreePy(data=as.dist(mdist), method="average")
  res <- tree_py(mdist, method="average")
  
  tree_filename <- file.path(
      data_dir,
      "gleason6.1/results/",
      "GL6.1smear1bpFisherTreePyP.txt")
  expect_true(file.exists(tree_filename))
  tree <- load_table(tree_filename)
  
  expect_equal(res[, 'index1'], tree[, 'index1'])
  expect_equal(res[, 'index2'], tree[, 'index2'])
  expect_equal(res[, 'height'], tree[, 'height'])
  
})