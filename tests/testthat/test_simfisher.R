context("check simfisher works as expected")

test_that("simfisher_wrapper works as expected", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  pinmat_filename <- file.path(data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  pinmat_df <- load_table(pinmat_filename)

  pins_filename <- file.path(data_dir, "gleason6.1/results/GL6.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  pins_df <- load_table(pins_filename)

  res <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=2, nsweep=200)

  expect_true(!is.null(res))

  true_res <- res$true

  true_filename <- file.path(data_dir, "gleason6.1/results/GL6.1trueP.txt")
  expect_true(file.exists(true_filename))
  true_mat <- as.matrix(scan(true_filename))

  expect_equal(true_mat, true_res, tolerance=1e-7)

  sim_res <- res$sim
  # print("sim_res>")
  # print(dim(sim_res))

  sim_filename <- file.path(data_dir, "gleason6.1/test_results/sim_fisher_test_2_of_500_by_200_simP.txt.gz")
  expect_true(file.exists(sim_filename))
  sim_mat <- as.matrix(scan(sim_filename))
  # print("sim_mat>")
  # print(dim(sim_mat))
  size <- ncol(sim_res)*nrow(sim_res)
  # print("size>")
  # print(size)

  sim_mat <- sim_mat[1:size, 1]
  # print("sim_mat>")
  # print(dim(sim_mat))
  dim(sim_res) <- dim(sim_mat)
  print(dim(sim_mat))
  print(dim(sim_res))

  sim_res <- sort(sim_res)
  sim_mat <- sort(sim_mat)

  expect_equal(sim_mat, sim_res, tolerance=1e-2)

})
