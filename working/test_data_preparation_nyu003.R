devtools::load_all()

data_dir <- Sys.getenv("SGAINS_DATA")
assertthat::assert_that(file.exists(data_dir))

nyu003_varbin_dir <- file.path(data_dir, "nyu003/varbin_orig")
nyu003_scgv_dir <- file.path(data_dir, "nyu003/scgv")
assertthat::assert_that(file.exists(nyu003_varbin_dir))
assertthat::assert_that(file.exists(nyu003_scgv_dir))

CJA3182_varbin_df <- load_table(file.path(nyu003_varbin_dir, "CJA3182.varbin.20k.txt.gz"))
assertthat::are_equal(nrow(CJA3182_varbin_df), 20000)

segment_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt"))
assertthat::are_equal(nrow(segment_df), 19943)

colnames(segment_df)

CJA3182_segment_df <- segment_df[c("chrom", "chrompos", "abspos", "CJA3182")]
colnames(CJA3182_segment_df)
dim(CJA3182_segment_df)

seg_filename <- file.path(data_dir, "nyu003/segmented/CJA3182.seg.txt")
write.table(CJA3182_segment_df, file=seg_filename, row.names = F, col.names = T, quote = F, sep="\t")

ratio_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.lowratio.quantal.R.ratio.txt"))
dim(ratio_df)
colnames(ratio_df)
assertthat::are_equal(nrow(ratio_df), 19943)

CJA3182_ratio_df <- ratio_df[c("chrom", "chrompos", "abspos", "CJA3182")]
dim(CJA3182_ratio_df)

ratio_filename <- file.path(data_dir, "nyu003/segmented/CJA3182.ratio.txt")
write.table(CJA3182_ratio_df, file=ratio_filename, row.names = F, col.names = T, quote = F, sep="\t")

gc_df <- load_table(file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz"))
assertthat::are_equal(nrow(gc_df), 20000)

bins_df <- augment_bins_sample(CJA3182_bins_df, gc_df)
res_df <- cbs.segment_1(bins_df, gc_df, alpha = 0.05, nperm = 1000, undo.SD = 1.0, min.width = 5)

head(res_df)
head(CJA3182_segment_df)
head(CJA3182_ratio_df)
m <- match(res_df$abspos, CJA3182_ratio_df$abspos, nomatch = F)
length(m)
df <- res_df["ratio.quantal",m]
head(res_df)
vals <- res_df$ratio.quantal[m]
length(vals)
abs(vals - CJA3182_ratio_df$CJA3182) < 1e-6

sum(abs(vals - CJA3182_ratio_df$CJA3182) < 1e-6)
