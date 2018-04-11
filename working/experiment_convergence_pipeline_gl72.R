# 
# Author: lubo
# Created: Mar 28, 2018
###############################################################################
library(devtools)
library(testthat)
library(futile.logger)

devtools::load_all()


get_test_sim_filename <- function(dirname, count, nsim, nsweep) {
  filename <- paste(
      "sim_fisher_test_",
      count, "_of_", nsim, "_by_", nsweep,
      "_simP.txt.gz",sep="")
  return(file.path(dirname, filename))
}


flog.threshold(DEBUG)
flog.debug("experimental convergence pipeline for GL7.2")

data_dir <- Sys.getenv("SGAINS_DATA")    

gc_filename <- file.path(
    data_dir, 
    "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
gc_df <- load_table(gc_filename)
gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)
gc_df <- gc_df[-badbins.20k$V1, ]
expect_equal(nrow(gc_df), 19943)

seg_filename <- file.path(
    data_dir,
    "gleason7.2/uber/uber.hg19.nyu007.GL7.2.20k.seg.quantal.R.txt.gz")
expect_true(file.exists(seg_filename))

segment_df <- load_table(seg_filename)
# segment_df <- filter_good_columns(segment_df, cells, skip=3)
# expect_equal(ncol(segment_df), length(cells)+3)
expect_equal(nrow(segment_df), 19943)


rat_filename <- file.path(
    data_dir,
    "gleason7.2/uber/uber.hg19.nyu007.GL7.2.20k.lowratio.quantal.R.txt.gz")
expect_true(file.exists(rat_filename))

ratio_df <- load_table(rat_filename)
# ratio_df <- filter_good_columns(ratio_df, cells, skip=3)
# expect_equal(ncol(ratio_df), length(cells)+3)
expect_equal(nrow(ratio_df), 19943)

res_dir <- file.path(data_dir, "gleason7.2/results")
true_filename <- file.path(
    res_dir, "nyu007.GL7.2trueP.txt")
pinmat_filename <- file.path(
    res_dir, "nyu007.GL7.2smear1bpPinMat.txt")
pins_filename <- file.path(
    res_dir, "nyu007.GL7.2smear1bpPins.txt")

expect_true(file.exists(true_filename))
expect_true(file.exists(pinmat_filename))
expect_true(file.exists(pins_filename))

true_pv <- scan(true_filename)
pinmat_df <- load_table(pinmat_filename)
pins_df <- load_table(pins_filename)

cell <- uber_cells(pinmat_df, skip=0)[,1]
cells <- cell

sim_dirname <- file.path(
    data_dir,
    "gleason7.2/out_many")
flog.debug("all data loaded...")


for(nsim in c(10)) {
  flog.debug("nsim=%s",nsim)

  sim_filename <- get_test_sim_filename(sim_dirname, nsim=500, 500, 200)
  sim_pv <- scan(sim_filename)
  sim_pv <- sim_pv[!is.na(sim_pv)]
  expect_false(any(is.na(sim_pv)))
  
  flog.debug("simulation loaded")
  
  mfdr <- fisher_fdr(true_pv, sim_pv, cells)
  mdist <- fisher_dist(true_pv, cells)
  
  hc <- hclust_tree(pinmat_df, mfdr, mdist)
  tree_df <- tree_py(mdist, method='average')
  hc <- find_clones(hc)
  subclones <- find_subclones(hc, pinmat_df, pins_df, nsim=nsim, lmmax=0.001)
  
  out_dir <- file.path(
      data_dir,
      "gleason7.2/conv_03", 
      paste("subGL_72_",nsim,sep=""))
  if(!file.exists(out_dir)) {
    dir.create(out_dir)
  }
  filenames <- case_filenames(out_dir, "sub_GL7_2")

  save_table(filenames$cells, data.frame(cell))
  save_table(filenames$seg, segment_df)
  save_table(filenames$ratio, ratio_df)
  save_table(filenames$featuremat, pinmat_df)
  save_table(filenames$features, pins_df)
  save_table(filenames$tree, tree_df)
  save_table(filenames$clone, subclones)
}
