context("chrom_numeric works properly")

test_that("chrom_numeric works for chr1 and chr2", {
  chrom <- c("chr1", "chr1", "chr2", "chr2")
  bin_mat <- data.frame(chrom)

  res <- chrom_numeric(bin_mat$chrom)
  head(res)
  expect_that(res, equals(c(1,1,2,2)))
})


test_that("chrom_numeric works for chrX and chrY", {
  chrom <- c("chrX", "chrX", "chrY", "chrY")
  bin_mat <- data.frame(chrom)

  res <- chrom_numeric(bin_mat$chrom)
  head(res)
  expect_that(res, equals(c(23,23,24,24)))
})
