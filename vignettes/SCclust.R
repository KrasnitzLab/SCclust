## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library("SCclust")

## ------------------------------------------------------------------------
library(futile.logger)
flog.threshold(ERROR)

## ---- results = "asis"---------------------------------------------------
gc_df <- read.csv("T10data/hg19_R50_B20k_bins_boundaries.txt", header = T, sep='\t')
gc_df$chrom.numeric <- chrom_numeric(gc_df$bin.chrom)
knitr::kable(head(gc_df))

## ---- results = "asis"---------------------------------------------------
dim(gc_df)

## ---- results = "asis"---------------------------------------------------
cytobands <- read.csv("T10data/cytoBand.txt", header = F, sep='\t')
knitr::kable(head(cytobands))

## ---- results = "asis"---------------------------------------------------
centroareas <- calc_centroareas(cytobands)
knitr::kable(head(centroareas))

## ---- results = "asis"---------------------------------------------------
sample_df <- read.csv("T10data/varbin/SRR052047.varbin.20k.txt", header=T, sep='\t')
knitr::kable(head(sample_df))

## ---- results = "asis"---------------------------------------------------
sample_df <- read.csv("T10data/varbin/SRR052148.varbin.20k.txt", header=T, sep='\t')
knitr::kable(head(sample_df))

## ---- results = "asis"---------------------------------------------------
centrobins <- calc_regions2bins(gc_df, centroareas)
length(centrobins)

## ---- results = "asis"---------------------------------------------------
gc_df <- gc_df[-centrobins, ]
dim(gc_df)

## ---- results = "asis"---------------------------------------------------
varbin_files <- varbin_input_files("T10data/varbin", "*.varbin.20k.txt")
knitr::kable(head(varbin_files))

## ------------------------------------------------------------------------
varbin_files <- varbin_files[seq(10),]
dim(varbin_files)

## ------------------------------------------------------------------------
res <- segment_varbin_files(varbin_files, gc_df, centrobins)

## ---- results = "asis"---------------------------------------------------
knitr::kable(head(res$seg))

## ---- results = "asis"---------------------------------------------------
knitr::kable(head(res$ratio))

## ------------------------------------------------------------------------
filenames <- case_filenames("T10data/results", "NavinT10")

## ------------------------------------------------------------------------
save_table(filenames$seg, res$seg)
save_table(filenames$ratio, res$ratio)

## ------------------------------------------------------------------------
cells <- uber_cells(res$seg)$cells
save_table(filenames$cells, data.frame(cell=cells))

## ------------------------------------------------------------------------
dir("T10data/results")

## ------------------------------------------------------------------------
pins <- calc_pinmat(gc_df, res$seg, dropareas=centroareas)
pinmat_df <- pins$pinmat
pins_df <- pins$pins
save_table(filenames$featuremat, pinmat_df)
save_table(filenames$features, pins_df)

## ------------------------------------------------------------------------
dir("T10data/results")

## ------------------------------------------------------------------------
fisher <- sim_fisher_wrapper(
      pinmat_df, pins_df, njobs=30, nsim=150, nsweep=10)
true_pv <- fisher$true
sim_pv <- fisher$sim

## ------------------------------------------------------------------------
# devtools::load_all()
# flog.threshold(DEBUG)
mfdr <- fisher_fdr(true_pv, sim_pv, cells)
mdist <- fisher_dist(true_pv, cells)

## ------------------------------------------------------------------------
save_mat(filenames$true_pv, true_pv)
save_mat(filenames$sim_pv, sim_pv)

## ------------------------------------------------------------------------
hc <- hclust_tree(pinmat_df, mfdr, mdist)
tree_df <- tree_py(mdist, method='average')

## ------------------------------------------------------------------------
save_table(filenames$tree, tree_df)

## ------------------------------------------------------------------------
hc <- find_clones(hc)
subclones <- find_subclones(hc, pinmat_df, pins_df, nsim=nsim)

## ------------------------------------------------------------------------
save_table(filenames$clone, subclones)

