context("check segmenting works as expected in nyu003 case")


test_that("check that CJA3182 is segmented properly", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  expect_true(file.exists(data_dir))

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  expect_equal(nrow(gc_df), 20000)

  nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin_orig")
  nyu003_scgv_dir <- file.path(data_dir, "nyu003/scgv")
  expect_true(file.exists(nyu003_varbin_dir))
  expect_true(file.exists(nyu003_scgv_dir))

  CJA3182_varbin_df <- load_table(file.path(nyu003_varbin_dir, "CJA3182.varbin.20k.txt.gz"))
  expect_equal(nrow(CJA3182_varbin_df), 20000)

  # segment_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt"))
  # expect_equal(nrow(segment_df), 19943)

  input_dir <- "/home/lubo/Work/sgains/data/bowtie-0.12.7/nyu003_R50_B20k_b0.12.7/varbin"
  suffix_pattern <- "\\.varbin.20k.txt$"

  files_list = varbin_input_files(input_dir, suffix_pattern)

  # print(head(files_list))
  # print(dim(files_list))

  expect_that(dim(files_list)[[1]], equals(317))
  expect_that(dim(files_list)[[2]], equals(3))

})
