
context("file utils work as expected")


test_that("we can load uncompressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)

  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 5)

})

test_that("we can load gzip compressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt.gz"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)

  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 5)

})

test_that("we can load bzip2 compressed tables", {
  input_dir <- "fixtures"
  filename <- "table.txt.bz2"

  filename <- file.path(input_dir, filename)
  expect_true(file.exists((filename)))

  df <- load_table(filename)

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

  expect_equal(filenames$cells, "output/proba.cells.txt")
  expect_equal(filenames$seg, "output/proba.seg.txt")
  expect_equal(filenames$ratio, "output/proba.ratio.txt")
  expect_equal(filenames$clone, "output/proba.clone.txt")

})
