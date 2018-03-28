# 
# Author: lubo
# Created: Mar 28, 2018
###############################################################################


test_that("experimental pipeline works as expected for GL6.1", {
  print(Sys.time())

  data_dir <- Sys.getenv("SGAINS_DATA")    
  cells_filename <- file.path(
      data_dir,
      "gleason6.1/goodcells23.txt")
  expect_true(file.exists(cells_filename))
  cells <- read.table(cells_filename, header=F, as.is=T)
  cells <- cells[, 1]
  print(cells)
  print(length(cells))
  out_dir <- "tmp"
  filenames <- case_filenames(out_dir, "sub_GL6_1")
  
  gc_filename <- file.path(
      data_dir, 
      "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
  gc_df <- load_table(gc_filename)
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)
  gc_df <- gc_df[-badbins.20k$V1, ]
  expect_equal(nrow(gc_df), 19943)

  print(Sys.time())
  seg_filename <- file.path(
      data_dir,
      "gleason6.1/uber/uber.hg19.GL6.1.20k.seg.quantal.R.txt.gz")
  expect_true(file.exists(seg_filename))
  
  segment_df <- load_table(seg_filename)
  segment_df <- filter_good_columns(segment_df, cells, skip=3)
  expect_equal(ncol(segment_df), length(cells)+3)
  expect_equal(nrow(segment_df), 19943)
  print("writing seg...")
  save_table(filenames$seg, segment_df)
  
  
  rat_filename <- file.path(
      data_dir,
      "gleason6.1/uber/uber.hg19.GL6.1.20k.lowratio.quantal.R.txt.gz")
  expect_true(file.exists(rat_filename))
  
  ratio_df <- load_table(rat_filename)
  ratio_df <- filter_good_columns(ratio_df, cells, skip=3)
  expect_equal(ncol(ratio_df), length(cells)+3)
  expect_equal(nrow(ratio_df), 19943)
  print("writing ratio...")
  save_table(filenames$ratio, ratio_df)
  print(Sys.time())

  cell <- uber_cells(segment_df, skip=3)[,1]
  save_table(filenames$cells, data.frame(cell))
  
  
  pins <- calc_pinmat(gc_df, segment_df, dropareas=dropareas)
  pinmat_df <- pins$pinmat
  pins_df <- pins$pins
  save_table(filenames$featuremat, pinmat_df)
  save_table(filenames$features, pins_df)
  
  print(Sys.time())
  fisher <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=5, nsweep=10)
  true_pv <- fisher$true
  sim_pv <- fisher$sim
  print(Sys.time())
  
  mfdr <- fisher_fdr(true_pv, sim_pv, cells)
  mdist <- fisher_dist(true_pv, cells)
  
  hc <- hclust_tree(pinmat_df, mfdr, mdist)
  
  tree_df <- tree_py(mdist, method='average')
  save_table(filenames$tree, tree_df)

  hc <- find_clones(hc)
  subclones <- find_subclones(hc, pinmat_df, pins_df)
  save_table(filenames$clone, subclones)

})