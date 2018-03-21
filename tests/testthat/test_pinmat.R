context("check pinmat works as expected")

test_that("calc_segments_short works as expected on CJA3182 and CJA3183", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1short20k.txt.gz")
  expect_true(file.exists(short_file))

  short_df <- load_table(short_file)

  cell_file <- file.path(data_dir, "nyu003/CJA3182_3.uber/CJA3182_3.seg.txt.gz")
  expect_true(file.exists(cell_file))
  segment_df <- load_table(cell_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]

  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_segments_short(augment_df, segment_df, homoloss=0)

  for(cell in c("CJA3182", "CJA3183")) {
    t1_df <- short_df[short_df$profid == cell,]
    t2_df <- result_df[result_df$profid == cell,]

    expect_equal(nrow(t1_df), nrow(t2_df))
    expect_true(all(t1_df$chrom == t2_df$chrom))
    expect_true(all(t1_df$chromstart == t2_df$chromstart))
    expect_true(all(t1_df$absstart == t2_df$absstart))
    expect_true(all(t1_df$chromend == t2_df$chromend))
    expect_true(all(t1_df$absend == t2_df$absend))
    expect_true(all(t1_df$segstart == t2_df$segstart))
    expect_true(all(t1_df$segend == t2_df$segend))
    expect_true(all(t1_df$segvals == t2_df$segvals))
    expect_true(all(t1_df$cvals == t2_df$cvals))
  }
})

test_that("calc_ploidies works as expected on CJA3182 and CJA3183", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  ploidies_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1ploidies.txt.gz")
  expect_true(file.exists(ploidies_file))

  ploidies_df <- load_table(ploidies_file)

  segment_file <- file.path(data_dir, "nyu003/CJA3182_3.uber/CJA3182_3.seg.txt.gz")
  expect_true(file.exists(segment_file))
  segment_df <- load_table(segment_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]
  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_ploidies(augment_df, segment_df)

  expect_equal(colnames(ploidies_df), colnames(result_df))
  expect_true(all(ploidies_df["CJA3182", ] == result_df["CJA3182", ]))
  expect_true(all(ploidies_df["CJA3183", ] == result_df["CJA3183", ]))

})


test_that("calc_segments_short works as expected on CJA3182, CJA3183, CJA3204, good_only=F", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1short20k.txt.gz")
  expect_true(file.exists(short_file))

  short_df <- load_table(short_file)

  segment_file <- file.path(data_dir, "nyu003/CJA3182_183_204.uber/CJA3182_183_204.seg.txt.gz")
  expect_true(file.exists(segment_file))
  segment_df <- load_table(segment_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]

  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_segments_short(augment_df, segment_df, homoloss=0)

  for(cell in c("CJA3182", "CJA3183", "CJA3204")) {
    t1_df <- short_df[short_df$profid == cell,]
    t2_df <- result_df[result_df$profid == cell,]

    expect_equal(nrow(t1_df), nrow(t2_df))
    expect_true(all(t1_df$chrom == t2_df$chrom))
    expect_true(all(t1_df$chromstart == t2_df$chromstart))
    expect_true(all(t1_df$absstart == t2_df$absstart))
    expect_true(all(t1_df$chromend == t2_df$chromend))
    expect_true(all(t1_df$absend == t2_df$absend))
    expect_true(all(t1_df$segstart == t2_df$segstart))
    expect_true(all(t1_df$segend == t2_df$segend))
    expect_true(all(t1_df$segvals == t2_df$segvals))
    expect_true(all(t1_df$cvals == t2_df$cvals))
  }
})


test_that("calc_segments_short works as expected on CJA3182, CJA3183, CJA3204, good_only=T", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1short20k.txt.gz")
  expect_true(file.exists(short_file))

  short_df <- load_table(short_file)

  segment_file <- file.path(data_dir, "nyu003/CJA3182_183_204.uber/CJA3182_183_204.seg.txt.gz")
  expect_true(file.exists(segment_file))
  segment_df <- load_table(segment_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]

  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_segments_short(augment_df, segment_df)

  expect_true(nrow(result_df[result_df$profid == "CJA3204", ]) == 0)

  for(cell in c("CJA3182", "CJA3183")) {
    t1_df <- short_df[short_df$profid == cell,]
    t2_df <- result_df[result_df$profid == cell,]

    expect_equal(nrow(t1_df), nrow(t2_df))
    expect_true(all(t1_df$chrom == t2_df$chrom))
    expect_true(all(t1_df$chromstart == t2_df$chromstart))
    expect_true(all(t1_df$absstart == t2_df$absstart))
    expect_true(all(t1_df$chromend == t2_df$chromend))
    expect_true(all(t1_df$absend == t2_df$absend))
    expect_true(all(t1_df$segstart == t2_df$segstart))
    expect_true(all(t1_df$segend == t2_df$segend))
    expect_true(all(t1_df$segvals == t2_df$segvals))
    expect_true(all(t1_df$cvals == t2_df$cvals))
  }
})

test_that("filter_segments_short works as expected on CJA3182, CJA3183, CJA3204, with CJA3204 as eviltwin",  {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1short20k.txt.gz")
  expect_true(file.exists(short_file))

  short_df <- load_table(short_file)

  segment_file <- file.path(data_dir, "nyu003/CJA3182_183_204.uber/CJA3182_183_204.seg.txt.gz")
  expect_true(file.exists(segment_file))
  segment_df <- load_table(segment_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]
  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_segments_short(augment_df, segment_df, homoloss=0)
  res_df <- filter_evil_short(result_df, eviltwins = c("CJA3204"))
  expect_true(!is.null(res_df))
  expect_false("CJA3204" %in% colnames(res_df))
})

test_that("filter_segments_short works as expected on CJA3182, CJA3183, CJA3204, with dropareas",  {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_file <- file.path(data_dir, "nyu003/results/nyu003.benign.1short20k.txt.gz")
  expect_true(file.exists(short_file))

  short_df <- load_table(short_file)

  segment_file <- file.path(data_dir, "nyu003/CJA3182_183_204.uber/CJA3182_183_204.seg.txt.gz")
  expect_true(file.exists(segment_file))
  segment_df <- load_table(segment_file)

  gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
  assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
  gc_df$bin.chrom <- chrom_numeric(gc_df$bin.chrom)

  gc_df <- gc_df[-badbins.20k$V1, ]
  augment_df <- augment_gc(gc_df, segment_df)

  result_df <- calc_segments_short(augment_df, segment_df, homoloss=0)
  res_df <- filter_dropareas_short(result_df, dropareas = dropareas)
  expect_true(!is.null(res_df))

})


test_that("calc_smear_breakpoints works as expected", {

  data_dir <- Sys.getenv("SGAINS_DATA")
  short_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1shortGood20k.txt")
  expect_true(file.exists(short_filename))

  short_df <- load_table(short_filename)

  censored_df <- calc_censored_index(short_df, dropareas)
  smear_df <- calc_smear_breakpoints(short_df, censored_df)

  expect_false(is.null(smear_df))

  res <- calc_pinmat(short_df, smear_df)
  expect_false(is.null(res))


  pinmat_df <- res$pinmat

  pinmat_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPinMat.txt")
  expect_true(file.exists(pinmat_filename))
  orig_df <- load_table(pinmat_filename)

  expect_equal(pinmat_df, orig_df)

  pins_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPins.txt")
  expect_true(file.exists(pins_filename))
  orig_df <- load_table(pins_filename)

  pins_df <- res$pins
  expect_equal(pins_df, orig_df)

})
