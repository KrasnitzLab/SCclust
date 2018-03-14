# TODO: Add comment
# 
# Author: lubo
###############################################################################


context("check hclust_tree works as expected")

test_that("check hclust_tree works as expected", {
  data_dir <- Sys.getenv("SGAINS_DATA")

  pinmat_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")
  dist_filename <- file.path(
    data_dir, "gleason6.1/results/GL6.1smear1bpLog10FisherP.txt")
  fdr_filename <- file.path(
    data_dir, "gleason6.1/results/GL6.1smear1bpLog10FisherFDR.txt")

  expect_true(file.exists(dist_filename))
  expect_true(file.exists(fdr_filename))
  expect_true(file.exists(pinmat_filename))

  dist <- load_table(dist_filename)
  fdr <- load_table(fdr_filename)

  pinmat <- load_table(pinmat_filename)
  cellnames <- uber_cells(pinmat, skip=0)[,1]


  hc <- hclust_tree(pinmat, fdr, dist)
  
  hc_filename <- file.path(
    data_dir,
    "gleason6.1/results/",
    "GL6.1smear1bpLog10FisherHCP.rda")
  expect_true(file.exists(hc_filename))
  attach(hc_filename)
      
  print(search())
  print(ls(pos=2))

  HCP <- GL6.1smear1bpLog10FisherHCP

  expect_false(is.null(HCP$meanfdr))
  expect_false(is.null(HCP$nodesize))
  expect_false(is.null(HCP$sharing))
  expect_false(is.null(HCP$complexity))
  
  expect_equal(hc$merge, HCP$merge, tolerance=1e-6)
  expect_equal(hc$meanfdr, HCP$meanfdr, tolerance=1e-6)
  expect_equal(hc$nodesize, HCP$nodesize, tolerance=1e-6)
  expect_equal(hc$sharing, HCP$sharing, tolerance=1e-6)
  expect_equal(hc$complexity, HCP$complexity, tolerance=1e-6)
  
  expect_false(is.null(HCP$labellist))
  expect_false(is.null(HCP$leaflist))
  
  expect_equal(hc$labellist, HCP$labellist)
  expect_equal(hc$leaflist, HCP$leaflist)

})

