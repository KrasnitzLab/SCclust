context("utils work as expected")

test_that("dropareas are constructed", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  cyto_file <- file.path(data_dir, "annot/cytoBandHG19.txt")
  expect_true(file.exists(cyto_file))

  cyto_data <- read.table(cyto_file, header=F,as.is=T)

  result <- calc_centroareas(cyto_data)
  expect_equal(dim(result)[1], 24)
  expect_true(all(result$chrom == dropareas$chrom))
  expect_true(all(result$from == dropareas$from))
  expect_true(all(result$to == dropareas$to))
})


test_that("badbins are constructed", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  cyto_file <- file.path(data_dir, "annot/cytoBandHG19.txt")
  expect_true(file.exists(cyto_file))

  cyto_data <- read.table(cyto_file, header=F,as.is=T)

  centroareas <- calc_centroareas(cyto_data)

  gc_file <- file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
  expect_true(file.exists(gc_file))

  gc <- load_table(gc_file)

  result <- calc_badbins(gc, centroareas)
  print(result)

})
