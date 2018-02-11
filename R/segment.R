
numeric_chrom <- function(chroms) {
  chrom.numeric <- substring(chroms, 4)
  chrom.numeric[which(chroms == "chrX")] <- "23"
  chrom.numeric[which(chroms == "chrY")] <- "24"
  chrom.numeric <- as.numeric(chrom.numeric)

  return(chrom.numeric)
}


augment_bins_sample <- function(bins_df, gc_df) {
  bins_df$chrom <- numeric_chrom(gc_df$bin.chrom)
  bins_df$gc.content <- gc_df$gc.content

  return(bins_df)
}


lowess.gc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}
