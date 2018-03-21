# 
# Author: lubo
# Created: Mar 21, 2018
###############################################################################


context("check simfisher convergence")

get_test_sim_filename <- function(dirname, count, nsim, nsweep) {
    filename <- paste(
      "sim_fisher_test_",
      count, "_of_", nsim, "_by_", nsweep,
      "_simP.txt.gz",sep="")
    return(file.path(dirname, filename))
}

#test_that("simfisher convergence experiments", {
#
#  sim_dirname <- file.path(
#      "/home/lubo/Work/SCclust/working/out_many/",
#      "Gleason6.1")
#  
#  sim_filename <- get_test_sim_filename(sim_dirname, 1, 500, 200)
#  expect_true(file.exists(sim_filename))
#
#  sim1 <- scan(sim_filename)
#  print(length(sim1))
#  sim1 <- sim1[!is.na(sim1)]
#  print(length(sim1))
#  l0 = length(sim1)
#
#  for(i in 2:2) {
#    sim_filename <- get_test_sim_filename(sim_dirname, i, 500, 200)
#    expect_true(file.exists(sim_filename))
#    sim2 <- scan(sim_filename)
#    print(length(sim2))
#    sim2 <- sim2[!is.na(sim2)]
#    print(length(sim2))
#    
#    expect_equal(i*l0, length(sim2))
#    l1 = length(sim1)
#
#    expect_equal(sim1, sim2[1:l1], tolerance=1e-6)
#    sim1 <- sim2
#  }
#})


test_that("simfisher convergence experiments full", {
  data_dir <- Sys.getenv("SGAINS_DATA")
  
  true_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1trueP.txt")
  pinmat_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")
  
  expect_true(file.exists(true_filename))
  expect_true(file.exists(pinmat_filename))
  
  true_mat <- scan(true_filename)
  pinmat <- load_table(pinmat_filename)
  cellnames <- uber_cells(pinmat, skip=0)[,1]

  sim_dirname <- file.path(
          "/home/lubo/Work/SCclust/working/out_many/",
          "Gleason6.1")
  
  sim_filename <- get_test_sim_filename(sim_dirname, 1, 500, 200)
  expect_true(file.exists(sim_filename))

  sim <- scan(sim_filename)
  sim <- sim[!is.na(sim)]

  mfdr1 <- fisher_fdr(true_mat, sim, cellnames)
  
  outfile <- file("gl6_convergence.txt", "w")
  writeLines("iteration,diff", con=outfile, sep='\n')
  for(i in 2:100) {
      sim_filename <- get_test_sim_filename(sim_dirname, i, 500, 200)
      sim <- scan(sim_filename)
      sim <- sim[!is.na(sim)]

      mfdr2 <- fisher_fdr(true_mat, sim, cellnames)
      
      diff <- max(abs(data.matrix(mfdr1) - data.matrix(mfdr2)))
      line <- paste(i, diff, sep=',')
      print(line)
      writeLines(line, con=outfile, sep='\n')

      mfdr1 <- mfdr1
  }
  close(outfile)
})
