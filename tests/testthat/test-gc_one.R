context("gc_one works properly")

test_that("gc_one for constant bincounts produces constant lowratio == 1", {
  gc.content <- c(rep(0.49,5),rep(0.51,5))
  gc <- data.frame(gc.content)

  bincount <- rep(10, 10)
  bin_mat <- data.frame(bincount)

  normalized <- gc_one(bin_mat, gc)
  expect_that(normalized$lowratio, equals(rep(1,10)))
})


test_that("gc_one for 2 diff bincounts produces constant lowratio == 1", {
  gc.content <- c(rep(0.49,5),rep(0.51,5))
  gc <- data.frame(gc.content)

  bincount <- c(rep(10, 5), rep(9,5))
  bin_mat <- data.frame(bincount)

  normalized <- gc_one(bin_mat, gc)
  expect_that(normalized$lowratio, equals(rep(1,10)))
})


test_that("gc_one for 2 diff bincounts inside constant gc range produce proportional lowratio", {
  gc.content <- c(rep(0.01,2),rep(0.99,2))
  gc <- data.frame(gc.content)

  for(k in c(2,3,4,5,6)) {
    bincount <- c(k*100,100,100,100)
    bin_mat <- data.frame(bincount)

    normalized <- gc_one(bin_mat, gc)
    lr = normalized$lowratio
    expect_equal(lr[1]/lr[2], k, 1e-2)
  }
})

test_that("gc_one for 2 diff bincounts inside constant gc range produce proportional lowratio", {
  gc.content <- c(rep(0.01,2),rep(0.99,2))
  gc <- data.frame(gc.content)

  for(k in c(2,3,4,5,6)) {
    bincount <- c(100,k*100,100,100)
    bin_mat <- data.frame(bincount)

    normalized <- gc_one(bin_mat, gc)
    lr = normalized$lowratio
    expect_equal(lr[2]/lr[1], k, 1e-2)
  }
})

test_that("gc_one for 2 diff bincounts inside constant gc range produce proportional lowratio", {
  gc.content <- c(rep(0.01,2),rep(0.99,2))
  gc <- data.frame(gc.content)

  for(k in c(2,3,4,5,6)) {
    bincount <- c(100,100,k*100,100)
    bin_mat <- data.frame(bincount)

    normalized <- gc_one(bin_mat, gc)
    lr = normalized$lowratio
    expect_equal(lr[3]/lr[4], k, 1e-2)
  }
})

test_that("gc_one for 2 diff bincounts inside constant gc range produce proportional lowratio", {
  gc.content <- c(rep(0.01,2),rep(0.99,2))
  gc <- data.frame(gc.content)

  for(k in c(2,3,4,5,6)) {
    bincount <- c(100,100,100,k*100)
    bin_mat <- data.frame(bincount)

    normalized <- gc_one(bin_mat, gc)
    lr = normalized$lowratio
    expect_equal(lr[4]/lr[3], k, 1e-2)
  }
})
