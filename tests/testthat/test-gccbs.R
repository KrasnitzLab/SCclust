context("gcCBS_one works properly")

## TODO: Rename context
## TODO: Add more tests

test_that("gc_one works", {
  gc.content <- c(rep(0.1,5),rep(0.9,5))
  gc <- data.frame(gc.content)

  bincount <- rep(10, 10)
  bin_mat <- data.frame(bincount)

  normalized <- gc_one(bin_mat, gc)
  expect_that(normalized$lowratio, equals(rep(1,10)))
})
