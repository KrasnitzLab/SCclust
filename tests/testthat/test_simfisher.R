context("check simfisher works as expected")

test_that("simfisher_wrapper works as expected", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  pinmat_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  pinmat_df <- load_table(pinmat_filename)

  pins_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  pins_df <- load_table(pins_filename)

  res <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=500, nsweep=200)
  # res <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=2, nsweep=200)

  expect_true(!is.null(res))

  true_res <- res$true

  true_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1trueP.txt")
  expect_true(file.exists(true_filename))
  true_mat <- as.matrix(scan(true_filename))

  expect_equal(true_mat, true_res, tolerance=1e-7)

  sim_res <- res$sim

  sim_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1simP.txt")
  # sim_filename <- file.path(data_dir, "simfisher_simple/sim_fisher_test_2_simP.txt")
  expect_true(file.exists(sim_filename))
  sim_mat <- as.matrix(scan(sim_filename))

  dim(sim_mat) <- dim(sim_res)
  expect_equal(sim_mat, sim_res, tolerance=1e-3)

})
