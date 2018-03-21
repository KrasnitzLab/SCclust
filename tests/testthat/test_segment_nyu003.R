context("check segmenting works as expected in nyu003 case")


#test_that("check that CJA3182 is segmented properly", {
#
#  data_dir <- Sys.getenv("SGAINS_DATA")
#  expect_true(file.exists(data_dir))
#
#  gc_df <- load_table(
#      file.path(data_dir, 
#          "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
#  expect_equal(nrow(gc_df), 20000)
#
#  nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin_orig")
#  nyu003_scgv_dir <- file.path(data_dir, "nyu003/scgv")
#  expect_true(file.exists(nyu003_varbin_dir))
#  expect_true(file.exists(nyu003_scgv_dir))
#
#  CJA3182_varbin_df <- load_table(file.path(nyu003_varbin_dir, "CJA3182.varbin.20k.txt.gz"))
#  expect_equal(nrow(CJA3182_varbin_df), 20000)
#
#  # segment_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt"))
#  # expect_equal(nrow(segment_df), 19943)
#
##  input_dir <- "/home/lubo/Work/sgains/data/bowtie-0.12.7/nyu003_R50_B20k_b0.12.7/varbin"
##  suffix_pattern <- "\\.varbin.20k.txt$"
##
##  files_list = varbin_input_files(input_dir, suffix_pattern)
##
##  # print(head(files_list))
##  # print(dim(files_list))
##
##  expect_that(dim(files_list)[[1]], equals(317))
##  expect_that(dim(files_list)[[2]], equals(3))
#
#})


test_that("check that CJA3182 is segmented properly", {
        
    data_dir <- Sys.getenv("SGAINS_DATA")
    expect_true(file.exists(data_dir))
  
    gc_df <- load_table(
            file.path(data_dir, 
                    "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
    expect_equal(nrow(gc_df), 20000)
  
    nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin_orig")
    expect_true(file.exists(nyu003_varbin_dir))
    
    bin_df <- load_table(file.path(nyu003_varbin_dir, "CJA3182.varbin.20k.txt.gz"))
    expect_equal(nrow(bin_df), 20000)
    
    bin_df$chrom.numeric <- chrom_numeric(bin_df$chrom)
    bin_df <- cbs_segment_ratio(gc_df, bin_df)
    res_df <- cbs_segment_varbin(bin_df)

    CJA3182_varbin_dir <- file.path(data_dir, "nyu003/CJA3182/processed.v2")
    data_filename <- file.path(
        CJA3182_varbin_dir, 
        "CJA3182.hg19.20k.k50.varbin.data.txt")
    expect_true(file.exists(data_filename))
    data_df <- load_table(data_filename)
    
    expect_equal(res_df$chrom, data_df$chrom)
    expect_equal(res_df$chrompos, data_df$chrompos)
    expect_equal(res_df$abspos, data_df$abspos)
    expect_equal(res_df$ratio.quantal, data_df$ratio.quantal, tolerance=1e-6)
    expect_equal(res_df$seg.quantal, data_df$seg.quantal, tolerance=1e-2)
    gc(res_df)
    gc(data_df)
})


test_that("check that we can segment CJA3182 varbin produced by sgains", {
      
      data_dir <- Sys.getenv("SGAINS_DATA")
      expect_true(file.exists(data_dir))
    
      gc_df <- load_table(
      file.path(data_dir, 
          "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
      expect_equal(nrow(gc_df), 20000)
    
      nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin")
      expect_true(file.exists(nyu003_varbin_dir))
      
      bin_df <- load_table(file.path(nyu003_varbin_dir, "CJA3182.varbin.20k.txt.gz"))
      expect_equal(nrow(bin_df), 20000)
      
      bin_df$chrom.numeric <- chrom_numeric(bin_df$chrom)
      bin_df <- cbs_segment_ratio(gc_df, bin_df)
      
      res_df <- cbs_segment_varbin(bin_df)
      
      CJA3182_varbin_dir <- file.path(data_dir, "nyu003/CJA3182/processed.v2")
      data_filename <- file.path(
              CJA3182_varbin_dir, 
              "CJA3182.hg19.20k.k50.varbin.data.txt")
      expect_true(file.exists(data_filename))
      data_df <- load_table(data_filename)
      
      expect_equal(res_df$chrom.numeric, data_df$chrom)
      expect_equal(res_df$chrompos, data_df$chrompos)
      expect_equal(res_df$abspos, data_df$abspos)
      diff <- sum(abs(res_df$ratio.quantal - data_df$ratio.quantal) > 1e-3)
      error <- (1.0 * diff) / nrow(res_df)
      print(error)
      expect_true(error < 0.03)
      expect_equal(res_df$ratio.quantal, data_df$ratio.quantal, tolerance=1e-2)
      
      diff <- sum(abs(res_df$seg.quantal - data_df$seg.quantal) > 1e-3)
      error <- (1.0 * diff) / nrow(res_df)
      print(error)
      # expect_true(error < 0.003)
      expect_equal(res_df$seg.quantal, data_df$seg.quantal, tolerance=1e-2)

      gc(res_df)
      gc(data_df)
      
})
