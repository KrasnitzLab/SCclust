context("check fisher_fdr works as expected")

test_that("fisher_fdr works as expected", {
  data_dir <- Sys.getenv("SGAINS_DATA")

  true_filename <- file.path(
          data_dir, "gleason6.1/results/GL6.1trueP.txt")
  sim_filename <- file.path(
          data_dir, "gleason6.1/results/GL6.1simP.txt")
  pinmat_filename <- file.path(
          data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")

  expect_true(file.exists(true_filename))
  expect_true(file.exists(sim_filename))
  expect_true(file.exists(pinmat_filename))

  true_mat <- scan(true_filename)
  sim_mat <- scan(sim_filename)

  pinmat <- load_table(pinmat_filename)

  cellnames <- uber_cells(pinmat, skip=0)


  res <- fisher_fdr(true_mat, sim_mat, cellnames[,1])

  mfdr <- res$mat_fdr
  mdist <- res$mat_dist

  mdist_filename <- file.path(
          data_dir,
          "gleason6.1/results/",
          "GL6.1smear1bpLog10FisherP.txt")
  expect_true(file.exists(mdist_filename))
  mdist_bench <- load_table(mdist_filename)
  
  expect_equal(
          data.matrix(mdist), 
          data.matrix(mdist_bench), tolerance=1e-6)
  
  mfdr_filename <- file.path(
          data_dir,
          "gleason6.1/results/",
          "GL6.1smear1bpLog10FisherFDR.txt")
  expect_true(file.exists(mfdr_filename))
  mfdr_bench <- load_table(mfdr_filename)

  expect_equal(
          data.matrix(mfdr), 
          data.matrix(mfdr_bench), tolerance=1e-6)
})
