context("file utils work as expected")

test_that("varbin input file are constructed", {

  input_dir <- "/home/lubo/Work/sgains/data/bowtie-0.12.7/nyu003_R50_B20k_b0.12.7/varbin"
  suffix_pattern <- "\\.varbin.20k.txt$"

  files_list = varbin_input_files(input_dir, suffix_pattern)

  # print(head(files_list))
  # print(dim(files_list))

  expect_that(dim(files_list)[[1]], equals(317))
  expect_that(dim(files_list)[[2]], equals(3))

})

test_that("we can load uncompressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)
  # print(dim(df))

  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 5)

})

test_that("we can load gzip compressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt.gz"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)
  # print(dim(df))

  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 5)

})

test_that("we can load bzip2 compressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt.bz2"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)
  # print(dim(df))

  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 5)

})

test_that("we can extract uber cells names", {
  input_dir <- "fixtures"
  filename <- "table.txt.bz2"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)

  cells <- uber_cells(df)
  expect_equal(nrow(cells), 2)
})


test_that("we can extract uber cells names without skipping columns", {
  input_dir <- "fixtures"
  filename <- "table.txt.bz2"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)

  cells <- uber_cells(df, skip=0)
  expect_equal(nrow(cells), 5)
})


test_that("case filenames are build properly", {
  casename <- "proba"
  filenames <- case_filenames("output", casename)
  # print(filenames)

  expect_equal(filenames$cells, "output/proba.cells.txt")
  expect_equal(filenames$seg, "output/proba.seg.txt")
  expect_equal(filenames$ratio, "output/proba.ratio.txt")
  expect_equal(filenames$clone, "output/proba.clone.txt")

})
