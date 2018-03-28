# 
# Author: lubo
# Created: Mar 28, 2018
###############################################################################


context("check pinmat pipeline works as expected")


test_that("pinmat full build works as expected for nyu003", {
      
  data_dir <- Sys.getenv("SGAINS_DATA")

  pinmat_filename <- file.path(
      data_dir, "nyu003/results/nyu003.benign.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  pinmat_df <- load_table(pinmat_filename)
  
  cells <-uber_cells(pinmat_df, skip=0)$cells
  print(length(cells))
  
  gc_filename <- file.path(
      data_dir, 
      "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
  expect_true(file.exists(gc_filename))
  gc_df <- load_table(gc_filename)
  expect_equal(nrow(gc_df), 20000)

  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)
  gc_df <- gc_df[-badbins.20k$V1, ]
  expect_equal(nrow(gc_df), 19943)
  
  seg_filename <- file.path(
      data_dir,
      "nyu003/scgv/uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt")
  expect_true(file.exists(seg_filename))
  
  segment_df <- load_table(seg_filename)
  expect_equal(nrow(segment_df), 19943)

  augment_df <- augment_gc(gc_df, segment_df)

  short_df <- calc_segments_short(augment_df, segment_df, homoloss=0)
  short_df <- filter_good_short(short_df, cells)
  
  censored_index <- calc_censored_index(short_df, dropareas)
  print(length(censored_index))
  
  smear_df <- calc_smear_breakpoints(short_df, censored_index)
  
  res <- calc_pinmat_short(short_df, smear_df)
    
  pinmat_res <- res$pinmat
  expect_equal(pinmat_df, pinmat_res)
  
  pins_filename <- file.path(
      data_dir, "nyu003/results/nyu003.benign.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  pins_df <- load_table(pins_filename)
  
  pins_res <- res$pins
  expect_equal(pins_df, pins_res)
})

test_that("pinmat full build works as expected for GL6.1", {
  
  data_dir <- Sys.getenv("SGAINS_DATA")

  pinmat_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  pinmat_df <- load_table(pinmat_filename)
  
  cells <-uber_cells(pinmat_df, skip=0)$cells
  print(length(cells))
  
  gc_filename <- file.path(
          data_dir, 
          "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
  expect_true(file.exists(gc_filename))
  gc_df <- load_table(gc_filename)
  expect_equal(nrow(gc_df), 20000)

  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)
  gc_df <- gc_df[-badbins.20k$V1, ]
  expect_equal(nrow(gc_df), 19943)
  
  seg_filename <- file.path(
          data_dir,
          "/gleason6.1/uber/uber.hg19.GL6.1.20k.seg.quantal.R.txt.gz")
  expect_true(file.exists(seg_filename))
  segment_df <- load_table(seg_filename)
  expect_equal(nrow(segment_df), 19943)

  augment_df <- augment_gc(gc_df, segment_df)

  short_df <- calc_segments_short(augment_df, segment_df, homoloss=0)
  short_df <- filter_good_short(short_df, cells)
  
  censored_index <- calc_censored_index(short_df, dropareas)
  print(length(censored_index))
  
  smear_df <- calc_smear_breakpoints(short_df, censored_index)
  
  res <- calc_pinmat_short(short_df, smear_df)
  
  pinmat_res <- res$pinmat  
  expect_equal(pinmat_df, pinmat_res)
  
  pins_filename <- file.path(
      data_dir, "gleason6.1/results/GL6.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  pins_df <- load_table(pins_filename)
  
  pins_res <- res$pins
  expect_equal(pins_df, pins_res)
  
})