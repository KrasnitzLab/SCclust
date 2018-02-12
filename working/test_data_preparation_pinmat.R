devtools::load_all()

data_dir <- Sys.getenv("SGAINS_DATA")
assertthat::assert_that(file.exists(data_dir))

nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin_orig")
nyu003_scgv_dir <- file.path(data_dir, "nyu003/scgv")
assertthat::assert_that(file.exists(nyu003_varbin_dir))
assertthat::assert_that(file.exists(nyu003_scgv_dir))

segment_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt.gz"))
assertthat::are_equal(nrow(segment_df), 19943)

gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
assertthat::assert_that(assertthat::are_equal(nrow(gc_df), 20000))
gc_df$bin.chrom <- numeric_chrom(gc_df$bin.chrom)

nobad_gc_df <- gc_df[-badbins.20k$V1,]
assertthat::assert_that(assertthat::are_equal(nrow(nobad_gc_df), 19943))

head(nobad_gc_df)
head(segment_df[,c(1:10)])

assertthat::assert_that(all(nobad_gc_df$bin.chrom == segment_df$chrom))
assertthat::assert_that(all(nobad_gc_df$bin.start == segment_df$chrompos))

