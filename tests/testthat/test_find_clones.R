# 
# Author: lubo
###############################################################################


context("check find_clones works as expected")

test_that("check find_clones works as expected on GL6.1 case", {
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
  hc <- hclust_tree(pinmat, fdr, dist)
  
  hc <- find_clones(hc)

  hcp_filename <- file.path(
        data_dir,
        "gleason6.1/results/",
        "GL6.1smear1bpLog10FisherHCP.rda")
  expect_true(file.exists(hcp_filename))
  attach(hcp_filename)
  HCP <- GL6.1smear1bpLog10FisherHCP

  expect_false(is.null(HCP$fdrthresh))
  expect_false(is.null(HCP$sharemin))
  expect_false(is.null(HCP$nshare))

  expect_false(is.null(HCP$clonenodes))
  expect_false(is.null(HCP$softclones))
  
  expect_equal(hc$fdrthresh, HCP$fdrthresh, tolerance=1e-6)
  
  # default sharemin changed from 0.9 to 0.85
  expect_equal(hc$sharemin, HCP$sharemin, tolerance=0.05001)
  expect_equal(hc$nshare, HCP$nshare, tolerance=1e-6)
  
  expect_equal(hc$clonenodes, HCP$clonenodes)
  expect_equal(hc$softclones, HCP$softclones)

})


test_that("check find_clones works as expected on GL7.2 case", {
        data_dir <- Sys.getenv("SGAINS_DATA")
        res_dir <- file.path(data_dir, "gleason7.2/results")
        pinmat_filename <- file.path(
          res_dir, "nyu007.GL7.2smear1bpPinMat.txt")
        dist_filename <- file.path(
                  res_dir, "nyu007.GL7.2smear1bpLog10FisherP.txt")
        fdr_filename <- file.path(
                  res_dir, "nyu007.GL7.2smear1bpLog10FisherFDR.txt")
      
        expect_true(file.exists(dist_filename))
        expect_true(file.exists(fdr_filename))
        expect_true(file.exists(pinmat_filename))
      
        dist <- load_table(dist_filename)
        fdr <- load_table(fdr_filename)
      
        pinmat <- load_table(pinmat_filename)
        hc <- hclust_tree(pinmat, fdr, dist)
        
        hc <- find_clones(hc)
      
        hcp_filename <- file.path(
                  res_dir,
                  "nyu007.GL7.2smear1bpLog10FisherHCP.rda")
        expect_true(file.exists(hcp_filename))
        attach(hcp_filename)
        HCP <- nyu007.GL7.2smear1bpLog10FisherHCP
      
        expect_false(is.null(HCP$fdrthresh))
        expect_false(is.null(HCP$sharemin))
        expect_false(is.null(HCP$nshare))
      
        expect_false(is.null(HCP$clonenodes))
        expect_false(is.null(HCP$softclones))
        
        expect_equal(hc$fdrthresh, HCP$fdrthresh, tolerance=1e-6)
        expect_equal(hc$sharemin, HCP$sharemin, tolerance=1e-6)
        expect_equal(hc$nshare, HCP$nshare, tolerance=1e-6)
        
        expect_equal(hc$clonenodes, HCP$clonenodes)
        expect_equal(hc$softclones, HCP$softclones)
      
    })
