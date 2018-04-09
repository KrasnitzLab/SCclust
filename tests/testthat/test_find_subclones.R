# 
# Author: lubo
# Created: Mar 14, 2018
###############################################################################


context("check find_subclones works as expected")

test_that("check find_subclones works as expected GL6.1", {
  data_dir <- Sys.getenv("SGAINS_DATA")
  
  hcp_filename <- file.path(
      data_dir,
      "gleason6.1/results/",
      "GL6.1smear1bpLog10FisherHCP.rda")
  expect_true(file.exists(hcp_filename))
  attach(hcp_filename)
  
  print(hcp_filename)
  print(substring(hcp_filename,first=1,last=nchar(hcp_filename)-4))
  
  HCP <- GL6.1smear1bpLog10FisherHCP
  # print(ls(HCP))
  
  expect_true(!is.null(HCP$softclones))
  expect_true(!is.null(HCP$leaflist))
  expect_true(!is.null(HCP$labellist))
  expect_true(!is.null(HCP$labels))
  
  
  pinmat_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  pinmat <- load_table(pinmat_filename)
  
  pins_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  pins <- load_table(pins_filename)
  
  subclones <- find_subclones(HCP, pinmat, pins)
  
  # print(subclones)
  # print(length(subclones))
  
  
  subclones_filename <- file.path(
      data_dir,
      "gleason6.1/results/",
      "GL6.1smear1bpFisherPcloneTracks.txt")
  expect_true(file.exists(subclones_filename))
  
  subclones_bench <- load_table(subclones_filename)
  
  # print(subclones)
  # print(subclones_bench)
  
  expect_equal(subclones$ID, subclones_bench$ID)
  expect_equal(subclones$clone, subclones_bench$clone)
  expect_equal(subclones$subclone, subclones_bench$subclone)
})


test_that("check find_subclones works as expected GL7.2", {
    data_dir <- Sys.getenv("SGAINS_DATA")
    res_dir <- file.path(data_dir, "gleason7.2/results")
    
    hcp_filename <- file.path(
        res_dir,
        "nyu007.GL7.2smear1bpLog10FisherHCP.rda")
    expect_true(file.exists(hcp_filename))
    
    attach(hcp_filename)
    
    print(hcp_filename)
    print(substring(hcp_filename,first=1,last=nchar(hcp_filename)-4))
    
    HCP <- nyu007.GL7.2smear1bpLog10FisherHCP
    # print(ls(HCP))
    
    expect_true(!is.null(HCP$softclones))
    expect_true(!is.null(HCP$leaflist))
    expect_true(!is.null(HCP$labellist))
    expect_true(!is.null(HCP$labels))
    
    
    pinmat_filename <- file.path(
        res_dir, "nyu007.GL7.2smear1bpPinMat.txt")
    expect_true(file.exists(pinmat_filename))
    pinmat <- load_table(pinmat_filename)
    
    pins_filename <- file.path(
        res_dir, "nyu007.GL7.2smear1bpPins.txt")
    expect_true(file.exists(pins_filename))
    pins <- load_table(pins_filename)
    
    subclones <- find_subclones(HCP, pinmat, pins)
    
    # print(subclones)
    # print(length(subclones))
    
    
    subclones_filename <- file.path(
            res_dir,
            "nyu007.GL7.2smear1bpFisherPcloneTracks.txt")
    expect_true(file.exists(subclones_filename))
    
    subclones_bench <- load_table(subclones_filename)
    
    # print(subclones)
    # print(subclones_bench)
    save_table("temp.txt", subclones)
    
    expect_equal(subclones$ID, subclones_bench$ID)
    expect_equal(subclones$clone, subclones_bench$clone)
    expect_equal(subclones$subclone, subclones_bench$subclone)
})
