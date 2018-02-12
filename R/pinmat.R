

augment_gc <- function(gc_df, df) {
  assertthat::assert_that(all(gc_df$bin.chrom == df$chrom))
  assertthat::assert_that(all(gc_df$bin.start == df$chrompos))

  augment_df <- cbind(gc_df[,c("bin.chrom", "bin.start", "bin.end")], df[,"abspos"])
  colnames(augment_df) <- c("chrom", "chromstart", "chromend", "absstart")

  augment_df <- cbind(
    augment_df,
    augment_df[,"absstart"] + (augment_df[,"chromend"] - augment_df[,"chromstart"])
  )
  colnames(augment_df) <- c("chrom", "chromstart", "chromend", "absstart", "absend")
  return(augment_df)
}

calc_segments_short <- function(gc_df, segment_df) {
  assertthat::assert_that(((typeof(gc_df$chrom) == "integer") | (typeof(gc_df$chrom) == "double")))
  assertthat::assert_that(typeof(segment_df$chrom) == "integer")
  assertthat::assert_that(all(colnames(gc_df) == c("chrom", "chromstart", "chromend", "absstart", "absend")))
  assertthat::assert_that(nrow(gc_df) == nrow(segment_df))

  a <- round(as.matrix(segment_df[,-(1:3)]))

  b <- a-rbind(matrix(nrow=1,ncol=ncol(a),data=0),a[-nrow(a),])

  bb <- (b!=0)
  bb[match(unique(segment_df[,"chrom"]),segment_df[,"chrom"]),] <- T

  segstarts <- row(a)[bb]
  segvals <- a[bb]
  profid <- dimnames(a)[[2]][cumsum(segstarts==1)]
  segends <- c(segstarts[-1]-1,nrow(a))
  segends[segends==0] <- nrow(a)
  segbins <- segends-segstarts+1

  segloc <- cbind(gc_df[segstarts,c("chrom","chromstart","absstart")],
                  gc_df[segends,c("chromend","absend")])
  tshort <- data.frame(I(profid), segloc, segstarts, segends, segbins, segvals)

  ploidies_df <- calc_ploidies(gc_df, segment_df, a)
  print(dim(tshort))
  print(dim(ploidies_df))

  tshort <- cbind(tshort, tshort[,"segvals"] - ploidies_df[tshort[,"profid"], "ploidychromod"])
  colnames(tshort)[ncol(tshort)] <- "cvals"

  return(tshort)
}

calc_ploidies <- function(gc_df, segment_df, a=NULL) {
  if(is.null(a)) {
    a <- round(as.matrix(segment_df[,-(1:3)]))
  }
  assertthat::assert_that(nrow(a) == nrow(segment_df))
  assertthat::assert_that(ncol(a) + 3 == ncol(segment_df))

  getmode<-function(x) as.numeric(names(which.max(tapply(X=x,INDEX=as.factor(x), FUN=length))))
  ploidymod<-apply(a[gc_df[,"chrom"]<23,],2,getmode)
  ploidymed<-apply(a[gc_df[,"chrom"]<23,],2,median)
  ploidychromod<-apply(
    a[gc_df[,"chrom"]<23,], 2, modeofmodes,
    otherlabel=gc_df[gc_df[,"chrom"]<23,"chrom"],
    tiebreaker=2,tiebreakerside="greater")
  homoloss<-colSums(!a[gc_df[,"chrom"]<23,])/
    sum(gc_df[,"chrom"]<23)
  ploidies <- cbind(ploidymed, ploidymod, ploidychromod, homoloss)
  return(ploidies)
}


