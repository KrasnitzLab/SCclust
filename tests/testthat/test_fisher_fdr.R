context("check fisher_fdr works as expected")

test_that("fisher_fdr works as expected", {
  data_dir <- Sys.getenv("SGAINS_DATA")

  true_filename <- file.path(data_dir, "gleason6.1/results/GL6.1trueP.txt")
  sim_filename <- file.path(data_dir, "gleason6.1/results/GL6.1simP.txt")
  pinmat_filename <- file.path(data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")

  expect_true(file.exists(true_filename))
  expect_true(file.exists(sim_filename))
  expect_true(file.exists(pinmat_filename))

  true_mat <- scan(true_filename)
  sim_mat <- scan(sim_filename)

  pinmat <- load_table(pinmat_filename)

  cellnames <- uber_cells(pinmat, skip=0)


  res <- fdr_fisherPV(true_mat, sim_mat, cellnames, graphic=FALSE)

  mfdr <- res$mat_fdr
  mdist <- res$mat_dist

  print(dim(mfdr))
  print(dim(mdist))

})
