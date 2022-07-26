context("divisive clustering works as expected")

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
    expect_equal(length(incidence), 2)

})


test_that("we can use squeeze_vector", {
    v <- c(0, 0, 1, 1, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(28))

    v <- c(1, 0)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(1))

    v <- c(0, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(2))

    v <- c(0, 0, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(4))

    v <- c(0, 0, 0, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(8))

    v <- c(0, 0, 0, 0, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(16))

    v <- c(1, 0, 0, 1)
    r <- squeeze_vector(v)
    expect_equal(r, as.raw(9))

})

test_that("we can build incidence table 3x3", {
    m <- matrix(c(
        0, 0, 1,
        0, 1, 0,
        1, 0, 0),
        nrow=3, ncol=3, byrow=T)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)
    expect_false(is.null(dim(i[[1]])))

    expect_true( all(as.raw(c(04, 02, 01)) == i[[1]]) )
    expect_true( all(as.raw(c(04, 02, 01)) == i[[2]]) )
})


test_that("we can build matrix from incidence table 3x3", {
    m <- matrix(c(
        0, 0, 1,
        0, 1, 0,
        1, 0, 0),
        nrow=3, ncol=3, byrow=T)

    i <- build_incidence_table(m)

    r <- incidence_table2matrix(i)

    expect_equal(ncol(r), 3)
    expect_equal(nrow(r), 3)

    expect_true( all(m == r))
})


test_that("we can build matrix from incidence table 13x9", {
    m <- matrix(c(
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1),
        ncol=9, byrow=T)

    i <- build_incidence_table(m)

    r <- incidence_table2matrix(i)

    expect_equal(ncol(r), 9)
    expect_equal(nrow(r), 13)

    expect_true( all(m == r))
})


test_that("we can build incidence table 3x4", {
    m <- matrix(c(
        0, 0, 0, 1,
        0, 0, 1, 0,
        0, 1, 0, 0),
        nrow=3, ncol=4, byrow=T)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)
    expect_false(is.null(dim(i[[1]])))

    expect_true( all(as.raw(c(00, 04, 02, 01)) == i[[1]]) )
    expect_true( all(as.raw(c(08, 04, 02)) == i[[2]]) )
})

test_that("we can use replicate incidence 3x4", {
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
    input_dir <- "fixtures"
    filename <- "diagonal9.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    m <- load_matrix(filename)

    flog.debug("rows: %s", nrow(m))
    flog.debug("cols: %s", ncol(m))
    i <- build_incidence_table(m)

    r <- replicate_incidence(i, 1)

    expect_true( all(i[[1]] == r[[1]]) )
    expect_true( all(i[[2]] == r[[2]]) )

})

test_that("we can use consolidate incidence 3x4", {

    m <- matrix(c(
        0, 0, 0, 1,
        0, 0, 1, 0,
        0, 1, 0, 0),
        nrow=3, ncol=4, byrow=T)

    i <- build_incidence_table(m)

    r <- consolidate_incidence(i, 1)

    expect_true( all(as.raw(c(04, 02, 01)) == r[[1]]) )
    expect_true( all(as.raw(c(04, 02, 01)) == r[[2]]) )

})


test_that("we can use consolidate incidence 5x5", {

    m <- matrix(c(
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        0, 0, 0, 1, 1,
        0, 0, 0, 1, 1),
        ncol=5, byrow=T)

    i <- build_incidence_table(m)

    r <- consolidate_incidence(i, 1)
    flog.debug("length(r)=%s", length((r)))

    expect_equal(length(r), 2)

    flog.debug(as.vector(r[[1]]))
    # flog.debug("r[[2]]= %s", r[[2]])

    expect_true( all(as.raw(c(0x7, 0x7, 0x7, 0x18, 0x18)) == r[[1]]) )
    expect_true( all(as.raw(c(0x7, 0x7, 0x7, 0x18, 0x18)) == r[[2]]) )

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


test_that("we can use contingencies 5x5", {

    m <- matrix(c(
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        0, 0, 0, 1, 1,
        0, 0, 0, 1, 1),
        ncol=5, byrow=T)

    i <- build_incidence_table(m)

    p <- squeeze_vector(c(T, T, T, F, F))
    r <- contingencies(i, p)
    flog.debug("result: %s", r)

    expect_true(all(names(r) == c("contables", "pmarginals")))

    ct <- r$contables
    expect_true(all(dim(ct) == c(5, 4)))

    expect_true(all(ct[1, ] == c(3, 0, 0, 2)))
    expect_true(all(ct[2, ] == c(3, 0, 0, 2)))
    expect_true(all(ct[3, ] == c(3, 0, 0, 2)))

    expect_true(all(ct[4, ] == c(0, 2, 3, 0)))
    expect_true(all(ct[5, ] == c(0, 2, 3, 0)))

    cm <- r$pmarginals
    flog.debug("pmarginals: %s, %s", cm[1], cm[2])

    expect_true(all(cm == c(3, 2)))

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


test_that("we can use mimax 5x5", {

    m <- matrix(c(
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        0, 0, 0, 1, 1,
        0, 0, 0, 1, 1),
        ncol=5, byrow=T)

    i <- build_incidence_table(m)

    flog.info("expected:  %s", showBits(as.raw(c(0x07))))
    flog.info("expected:  %s", showBits(as.raw(c(0x18))))

    set.seed(1)
    for(index in 1:8) {
        r <- mimax(i)
        flog.info("mi=%s; partition: %s", r$mi, showBits(r$partition))

        expect_true(
            all(as.raw(c(0x07)) == r$partition) || all(as.raw(c(0x18)) == r$partition))
    }

})


test_that("we can use mimax 13x9", {
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

    flog.info("expected:  %s", showBits(as.raw(c(0xe0, 0x01))))

    set.seed(1)

    for(index in 1:8) {
        r <- mimax(i)
        flog.info("mi=%s; partition: %s", r$mi, showBits(r$partition))
        expect_true(
            all(as.raw(c(0x1f, 0x00)) == r$partition) ||
            all(as.raw(c(0xe0, 0x01)) == r$partition) )
    }

})

test_that("we can use mpshuffle", {

    set.seed(1)

    m <- matrix(c(
        1, 1, 0,
        1, 1, 0,
        0, 0, 1),
        nrow=3, ncol=3, byrow=T)

    i <- build_incidence_table(m)

    expect_true( all(as.raw(c(0x3, 0x3, 0x4)) == i[[1]]) )
    expect_true( all(as.raw(c(0x3, 0x3, 0x4)) == i[[2]]) )

    dummyarg <- 1
    r <- mpshuffle(dummyarg, i, 1)

    expect_true( all(as.raw(c(0x3, 0x5, 0x2)) == r[[1]]) )
    expect_true( all(as.raw(c(0x3, 0x5, 0x2)) == r[[2]]) )

    rr <- mpshuffle(dummyarg, r, 1)

    expect_true( all(as.raw(c(0x3, 0x6, 0x1)) == rr[[1]]) )
    expect_true( all(as.raw(c(0x5, 0x3, 0x2)) == rr[[2]]) )

})


test_that("we can use mpshuffle 13x9", {

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
    set.seed(1)

    dummyarg <- 1
    r <- mpshuffle(dummyarg, i, 1)

    expect_true( all(as.raw(
        c(0x6f, 0x6f, 0x6f, 0x3f, 0x6f, 0x90, 0x90, 0x90, 0x90)) == r[[1]][1,]) )
    expect_true( all(as.raw(
        c(0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0x1f, 0x1f, 0x1f)) == r[[1]][2,]) )

})



test_that("we can use randomimax 13x9", {

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
    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)

    set.seed(1)
    r <- randomimax(i, swappars=test_swappars, maxempv=0.25, miobserved=11)
    # e <- c(-22.62314, -28.76079, -32.12585, -42.40447, -46.90315)
    e <- c(-22.62314, -27.37449, -33.51214, -38.01082, -42.76217)
    flog.info("expected: %s", print(e))
    flog.info("actual  : %s", print(r))

    expect_true(all(abs(r - e) <= 10e-3))

})


test_that("we can use mimain on clean 13x9", {

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

    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)

    r <- mimain(i,
        maxgens=2, maxempv=0.25,
        saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 9)

})


test_that("we can use mimain 13x9 again", {

    m <- matrix(c(
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1),
        ncol=9, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(i, maxempv=0.25, saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 9)

})



test_that("we can use mimain 26x9", {

    set.seed(1)
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
        0, 0, 0, 0, 0, 1, 1, 1, 1,

        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1

        ),
        ncol=9, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)
    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(i, maxempv=0.25, saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 9)

})


test_that("we can use mimain 5x5", {

    m <- matrix(c(
        1, 1, 1, 0, 0,
        1, 1, 1, 0, 0,
        1, 0, 0, 0, 0,
        0, 0, 0, 0, 1,
        0, 0, 0, 1, 1),
        ncol=5, byrow=T)

    i <- build_incidence_table(m)

    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)
    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(i, maxempv=0.25, saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 5)
    hc <- mimosa2hc(r)
    plot(hc)

})


test_that("we can use mimain on simpleIncidenceTable", {

    input_dir <- "fixtures"
    filename <- "simpleIncidenceTable.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    m <- load_matrix(filename)

    i <- build_incidence_table(m)

    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(
        i, maxgens=2, maxempv=0.25, saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 35)

})


test_that("we can use mimos2hc 13x9", {

    m <- matrix(c(
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 1, 1, 1, 1),
        ncol=9, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=10,
        burnin=50,
        niter=50,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(
        i, maxgen=5, maxempv=0.25,
        saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 9)
    expect_equal(class(r), "mimosa")

    hc<-mimosa2hc(r)
    expect_equal(class(r), "mimosa")

    plot(hc)
    expect_equal(class(hc), "hclust")

})

test_that("we can use mimain on simpleIncidenceTable", {

    input_dir <- "fixtures"
    filename <- "simpleIncidenceTable.txt"

    filename <- file.path(input_dir, filename)
    expect_true(file.exists((filename)))

    m <- load_matrix(filename)

    i <- build_incidence_table(m)

    test_swappars <- list(
        configs=10,
        burnin=100,
        niter=100,
        permeas=100,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(
        i, maxgens=3, maxempv=0.1, saspars=test_saspars, swappars=test_swappars)

    expect_equal(length(r$pathcode), 35)
    hc<-mimosa2hc(r)
    plot(hc)

})


test_that("we can use mimain on navin T10", {

    flog.debug("testing on navin T10 started")
    input_dir <- "fixtures"
    filename <- "hg19_navin_T10.featuremat.txt"

    filename <- file.path(input_dir, filename)
    # expect_true(file.exists((filename)))

    m <- load_table(filename)

    i <- build_incidence_table(m)

    test_swappars <- list(
        configs=50,
        burnin=50,
        niter=50,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0.0265,
        restarts=10,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
        suddenfreeze=F,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    # for(maxgens in c(3, 5, 7)) {
    set.seed(1)
    maxgens <- 7

    r <- mimain(
        i, maxgens=maxgens, maxempv=0.25,
        saspars=test_saspars,
        swappars=test_swappars,
        useCores=20)

    # expect_equal(length(r$pathcode), 95)

    hc<-mimosa2hc(r)
    plot(hc, main=paste("we can use mimosa2hc on navin T10; maxgens=", maxgens))
    # }

    # hcd <- as.dendrogram(hc)
    # plot(
    #     hcd,
    #     main="we can use mimain on navin T10",
    #     horiz=T)

})



test_that("we can use mimosa2hc 13x15", {

    m <- matrix(c(
        1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
        ncol=15, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=50,
        burnin=10,
        niter=10,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    maxgen <- 5
    set.seed(1)
    r <- mimain(
        i,
        maxgen=maxgen, maxempv=0.05,
        saspars=test_saspars,
        swappars=test_swappars, useCores=20)


    hc<-mimosa2hc(r)

    plot(hc, main=paste("we can use mimosa2hc 13x15; maxgen=", maxgen))
    expect_equal(class(r), "mimosa")

    expect_equal(length(r$pathcode), 15)
    expect_equal(class(r), "mimosa")

    expect_equal(class(hc), "hclust")

})


test_that("we can use mimain 4x4", {

    m <- matrix(c(
        1, 1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 1, 1),
        ncol=4, byrow=T)

    i <- build_incidence_table(m)

    test_swappars <- list(
        configs=5,
        burnin=5,
        niter=5,
        permeas=50,
        choosemargin=0.5)
    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)
    r <- mimain(
        i, maxempv=0.25,
        saspars=test_saspars,
        swappars=test_swappars, useCores=20)

    expect_equal(length(r$pathcode), 4)
    hc <- mimosa2hc(r)
    plot(hc, main="we can use mimain 4x4")

})


test_that("we can use mimos2hc 13x8", {

    m <- matrix(c(
        1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1),
        ncol=8, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=10,
        burnin=10,
        niter=10,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=F,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)

    Rprof()
    r <- mimain(
        i, maxgen=7, maxempv=0.25,
        saspars=test_saspars, swappars=test_swappars)
    Rprof(NULL)
    summaryRprof()

    expect_equal(length(r$pathcode), 8)
    expect_equal(class(r), "mimosa")

    hc<-mimosa2hc(r)
    expect_equal(class(r), "mimosa")

    plot(hc, main=paste("we can use mimos2hc 13x8; maxgen=7"))
    expect_equal(class(hc), "hclust")

})

test_that("we can use mimos2hc 13x8; suddenfreeze=T", {

    m <- matrix(c(
        1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 0, 0, 1, 1, 1, 1),
        ncol=8, byrow=T)

    i <- build_incidence_table(m)
    test_swappars <- list(
        configs=50,
        burnin=10,
        niter=10,
        permeas=50,
        choosemargin=0.5)

    test_saspars <- list(
    	minfreq=0,
        restarts=10,
    	suddenfreeze=T,
        cooler=1.12,
        acceptance=0.234,
        sweepspercycle=10,
        maxcycles=20,
        stopatfreezeout=T,
    	stopifnoprogress=2,
        epsilon=0.0001,
        betafudge=1.0)

    set.seed(1)

    Rprof()
    r <- mimain(
        i, maxgen=7, maxempv=0.25,
        saspars=test_saspars, swappars=test_swappars)
    Rprof(NULL)
    summaryRprof()

    expect_equal(length(r$pathcode), 8)
    expect_equal(class(r), "mimosa")

    hc<-mimosa2hc(r)
    expect_equal(class(r), "mimosa")

    plot(hc, main=paste("we can use mimos2hc 13x8; maxgen=7"))
    expect_equal(class(hc), "hclust")

})

# test_that("we can use GDP for 13x8 example", {

#     m <- matrix(c(
#         1, 1, 1, 1, 0, 0, 0, 0,
#         1, 1, 1, 1, 0, 0, 0, 0,
#         1, 1, 0, 0, 0, 0, 0, 0,
#         1, 1, 0, 0, 0, 0, 0, 0,
#         1, 1, 0, 0, 0, 0, 0, 0,
#         1, 0, 0, 0, 0, 0, 0, 1,
#         1, 0, 0, 0, 0, 0, 0, 1,
#         0, 0, 0, 0, 0, 0, 1, 1,
#         0, 0, 0, 0, 0, 0, 1, 1,
#         0, 0, 0, 0, 0, 0, 1, 1,
#         0, 0, 0, 0, 1, 1, 1, 1,
#         0, 0, 0, 0, 1, 1, 1, 1,
#         0, 0, 0, 0, 1, 1, 1, 1),
#         ncol=8, byrow=T)

#     i <- build_incidence_table(m)
#     test_swappars <- list(
#         configs=5000,
#         burnin=500,
#         permeas=500,
#         choosemargin=0.5)

#     test_saspars <- list(
#         restarts=10,
#         cooler=1.12,
#         acceptance=0.234,
#         sweepspercycle=10,
#         maxcycles=20,
#         stopatfreezeout=T,
#     	stopifnoprogress=2,
#         epsilon=0.0001,
#         betafudge=1.0)

#     for(maxgen in c(3, 4, 5, 6, 7)) {
#         set.seed(1)
#         r <- mimain(
#             i, maxgen=7, maxempv=0.25,
#             saspars=test_saspars, swappars=test_swappars)

#         expect_equal(length(r$pathcode), 8)
#         expect_equal(class(r), "mimosa")

#         hc<-mimosa2hc(r)
#         expect_equal(class(r), "mimosa")

#         plot(hc, main=paste("we can use GDP for 13x8 example; maxgen=", maxgen))
#         expect_equal(class(hc), "hclust")
#     }

# })
