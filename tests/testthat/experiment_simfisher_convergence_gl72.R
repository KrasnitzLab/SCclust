# 
# Author: lubo
# Created: Mar 21, 2018
###############################################################################


context("check simfisher convergence for GL7.2")

get_test_sim_filename <- function(dirname, count, nsim, nsweep) {
    filename <- paste(
      "sim_fisher_test_",
      count, "_of_", nsim, "_by_", nsweep,
      "_simP.txt.gz",sep="")
    return(file.path(dirname, filename))
}


test_that("simfisher convergence experiments full on GL7.2", {
  data_dir <- Sys.getenv("SGAINS_DATA")
  res_dir <- file.path(data_dir, "gleason7.2/results")
  true_filename <- file.path(
      res_dir, "nyu007.GL7.2trueP.txt")
  pinmat_filename <- file.path(
      res_dir, "nyu007.GL7.2smear1bpPinMat.txt")
  
  expect_true(file.exists(true_filename))
  expect_true(file.exists(pinmat_filename))
  
  true_mat <- scan(true_filename)
  pinmat <- load_table(pinmat_filename)
  cellnames <- uber_cells(pinmat, skip=0)[,1]

  sim_dirname <- file.path(
      data_dir,
      "gleason7.2/out_many")
  
  sim_filename <- get_test_sim_filename(sim_dirname, 1, 500, 200)
  expect_true(file.exists(sim_filename))

  sim <- scan(sim_filename)
  sim <- sim[!is.na(sim)]

  mfdr1 <- fisher_fdr(true_mat, sim, cellnames)
  
  outfile <- file("gl72_convergence.txt", "w")
  writeLines("iteration,diff,diff2", con=outfile, sep='\n')
  for(i in 2:500) {
      sim_filename <- get_test_sim_filename(sim_dirname, i, 500, 200)
      sim <- scan(sim_filename)
      sim <- sim[!is.na(sim)]

      mfdr2 <- fisher_fdr(true_mat, sim, cellnames)
      
      diff <- max(abs(data.matrix(mfdr1) - data.matrix(mfdr2)))
      diff2 <- data.matrix(mfdr1) - data.matrix(mfdr2)
      diff2 <- diff2 * diff2
      diff2 <- sqrt(sum(diff2))
      
      line <- paste(i, diff, diff2, sep=',')
      print(line)
      writeLines(line, con=outfile, sep='\n')

      mfdr1 <- mfdr2
  }
  close(outfile)
})
