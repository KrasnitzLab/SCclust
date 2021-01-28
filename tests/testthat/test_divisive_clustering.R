context("bit utils work as expected")

test_that("we can use setBit function", {
    flog.debug("setBits called: %s", setBits(1))

    expect_equal(1, setBits(1))
    expect_equal(1, setBits(2))
    expect_equal(2, setBits(3))
})


test_that("we can use load_mat", {

    input_dir <- "fixtures"
    filename <- "simpleIncidenceTable.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    m <- load_matrix(filename)

    flog.debug(class(m))
    flog.debug(dim(m))
    flog.debug(nrow(m))
    flog.debug(ncol(m))

    expect_equal(class(m), "matrix")
    expect_equal(nrow(m), 40)
    expect_equal(ncol(m), 35)
})


test_that("we can use load_incidence_table", {

    input_dir <- "fixtures"
    filename <- "simpleIncidenceTable.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    incidence <- load_incidence_table(filename)

    expect_equal(class(incidence), "incidencetable")
    flog.debug(length(incidence))
    expect_equal(length(incidence), 3)

})


test_that("we can use squeeze_vector", {
    v <- c(0, 0, 1, 1, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(28))

    v <- c(1, 0)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(1))

    v <- c(0, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(2))

    v <- c(0, 0, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(4))

    v <- c(0, 0, 0, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(8))

    v <- c(0, 0, 0, 0, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(16))

    v <- c(1, 0, 0, 1)
    r <- squeeze_vector(v)
    flog.debug(r)
    expect_equal(r, as.raw(9))

})

test_that("we can build incidence table 3x3", {
    # m <- matrix(c(
    #     0, 0, 0, 1, 1, 
    #     0, 0, 0, 1, 1,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0),
    #     nrow=5, ncol=5, byrow=T)

    m <- matrix(c(
        0, 0, 1, 
        0, 1, 0,
        1, 0, 0),
        nrow=3, ncol=3, byrow=T)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)
    flog.debug(i)
    expect_false(is.null(dim(i[[1]])))

    expect_true( all(as.raw(c(04, 02, 01)) == i[[1]]) )
    expect_true( all(as.raw(c(04, 02, 01)) == i[[2]]) )
})


test_that("we can build incidence table 3x4", {
    # m <- matrix(c(
    #     0, 0, 0, 1, 1, 
    #     0, 0, 0, 1, 1,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0),
    #     nrow=5, ncol=5, byrow=T)

    m <- matrix(c(
        0, 0, 0, 1, 
        0, 0, 1, 0,
        0, 1, 0, 0),
        nrow=3, ncol=4, byrow=T)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)
    flog.debug(i)
    expect_false(is.null(dim(i[[1]])))

    expect_true( all(as.raw(c(00, 04, 02, 01)) == i[[1]]) )
    expect_true( all(as.raw(c(08, 04, 02)) == i[[2]]) )
})

test_that("we can use replicate incidence 3x4", {
    # m <- matrix(c(
    #     0, 0, 0, 1, 1, 
    #     0, 0, 0, 1, 1,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0),
    #     nrow=5, ncol=5, byrow=T)

    m <- matrix(c(
        0, 0, 0, 1, 
        0, 0, 1, 0,
        0, 1, 0, 0),
        nrow=3, ncol=4, byrow=T)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)

    r <- replicate_incidence(i, 1)
    expect_true( all(as.raw(c(00, 04, 02, 01)) == r[[1]]) )
    expect_true( all(as.raw(c(08, 04, 02)) == r[[2]]) )

    r <- replicate_incidence(i, 2)
    expect_true( all(as.raw(c(00, 04, 02, 01)) == r[[1]]) )
    expect_true( all(as.raw(c(08, 04, 02)) == r[[2]]) )
})

test_that("we can use replicate incidence diagonal9", {
    # m <- matrix(c(
    #     0, 0, 0, 1, 1, 
    #     0, 0, 0, 1, 1,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0),
    #     nrow=5, ncol=5, byrow=T)

    # m <- matrix(c(
    #     0, 0, 0, 1, 
    #     0, 0, 1, 0,
    #     0, 1, 0, 0),
    #     nrow=3, ncol=4, byrow=T)

    input_dir <- "fixtures"
    filename <- "diagonal9.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    m <- load_matrix(filename)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)
    flog.debug(i)

    r <- replicate_incidence(i, 1)
    flog.debug(r[[1]])
    flog.debug(r[[2]])

    expect_true( all(i[[1]] == r[[1]]) )
    expect_true( all(i[[2]] == r[[2]]) )

})

test_that("we can use consolidate incidence 3x4", {
    # m <- matrix(c(
    #     0, 0, 0, 1, 1, 
    #     0, 0, 0, 1, 1,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0,
    #     1, 1, 1, 0, 0),
    #     nrow=5, ncol=5, byrow=T)

    m <- matrix(c(
        0, 0, 0, 1, 
        0, 0, 1, 0,
        0, 1, 0, 0),
        nrow=3, ncol=4, byrow=T)

    i <- build_incidence_table(m)

    flog.debug(i)

    r <- consolidate_incidence(i, 1)
    expect_true( all(as.raw(c(04, 02, 01)) == r[[1]]) )
    expect_true( all(as.raw(c(04, 02, 01)) == r[[2]]) )

    # r <- replicate_incidence(i, 2)
    # expect_true( all(as.raw(c(00, 04, 02, 01)) == r[[1]]) )
    # expect_true( all(as.raw(c(08, 04, 02)) == r[[2]]) )
})


test_that("we can use contingencies 3x3", {

    m <- matrix(c(
        0, 0, 1, 
        0, 1, 0,
        1, 0, 0),
        nrow=3, ncol=3, byrow=T)

    i <- build_incidence_table(m)

    p <- squeeze_vector(c(T, F, F))
    r <- contingencies(i, p)
    flog.debug("result: %s", r)

    expect_true(all(names(r) == c("contables", "pmarginals")))

    ct <- r$contables
    expect_true(all(dim(ct) == c(3, 4)))

    expect_true(all(ct[1, ] == c(0, 1, 1, 1)))
    expect_true(all(ct[2, ] == c(0, 1, 1, 1)))
    expect_true(all(ct[3, ] == c(1, 0, 0, 2)))

    cm <- r$pmarginals
    flog.debug("pmarginals: %s, %s", cm[1], cm[2])

    expect_true(all(cm == c(1, 2)))

})


test_that("we can use misum", {

    m <- matrix(c(
        1, 1, 0, 
        1, 1, 0,
        0, 0, 1),
        nrow=3, ncol=3, byrow=T)

    i <- build_incidence_table(m)

    p <- squeeze_vector(c(T, T, F))
    c <- contingencies(i, p)
    flog.debug("contingencies: %s", c)

    # expect_true(all(names(c) == c("contables", "pmarginals")))

    mi <- misum(c)
    flog.debug("result: %s", mi)
    expect_equal(mi, 0)
})


test_that("we can use miupdate", {

    m <- matrix(c(
        1, 1, 0, 
        1, 1, 0,
        0, 0, 1),
        nrow=3, ncol=3, byrow=T)

    i <- build_incidence_table(m)

    p <- squeeze_vector(c(T, T, F))
    c <- contingencies(i, p)
    flog.debug("contingencies: %s", c)

    ct <- c$contables
    expect_true(all(dim(ct) == c(3, 4)))

    expect_true(all(ct[1, ] == c(2, 0, 0, 1)))
    expect_true(all(ct[2, ] == c(2, 0, 0, 1)))
    expect_true(all(ct[3, ] == c(0, 1, 2, 0)))

    # cm <- r$pmarginals
    # flog.debug("pmarginals: %s, %s", cm[1], cm[2])


    # expect_true(all(names(c) == c("contables", "pmarginals")))

    r <- miupdate(1, i, c, p)

    ct <- c$contables
    expect_true(all(dim(ct) == c(3, 4)))

    expect_true(all(ct[1, ] == c(2, 0, 0, 1)))
    expect_true(all(ct[2, ] == c(2, 0, 0, 1)))
    expect_true(all(ct[3, ] == c(0, 1, 2, 0)))

})


test_that("we can use mimax", {

    m <- matrix(c(
        1, 1, 1, 1, 1, 0, 0, 0, 0, 
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1),
        ncol=9, byrow=T)

    i <- build_incidence_table(m)

    for(k in 1:15) {
        r <- mimax(i)
        showBits(r$partition)
        flog.debug("partition: %s", showBits(r$partition))
    }
})

