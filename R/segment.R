

augment_bins_sample <- function(bins_df, gc_df) {
  bins_df$chrom <- numeric_chrom(gc_df$bin.chrom)
  bins_df$gc.content <- gc_df$gc.content

  return(bins_df)
}


