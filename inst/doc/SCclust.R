## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library("SCclust")

## ---- results = "asis"---------------------------------------------------
gc_df <- read.csv("T10data/hg19_R50_B20k_bins_boundaries.txt", header = T, sep='\t')
knitr::kable(head(gc_df))

## ---- results = "asis"---------------------------------------------------
cytobands <- read.csv("T10data/cytoBand.txt", header = F, sep='\t')
knitr::kable(head(cytobands))

## ---- results = "asis"---------------------------------------------------
centroareas <- calc_centroareas(cytobands)
knitr::kable(head(centroareas, 5))

## ---- results = "asis"---------------------------------------------------
sample_df <- read.csv("T10data/varbin/SRR052047.varbin.20k.txt", header=T, sep='\t')
knitr::kable(head(sample_df))

## ---- results = "asis"---------------------------------------------------
sample_df <- read.csv("T10data/varbin/SRR052148.varbin.20k.txt", header=T, sep='\t')
knitr::kable(head(sample_df))

## ------------------------------------------------------------------------
# centrobins <- calc_regions2bins(gc_df, centroareas)

