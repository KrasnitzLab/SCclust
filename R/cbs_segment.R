# 
# Author: lubo
# Created: Mar 21, 2018
###############################################################################


lowess_gc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}


cbs_segment_ratio <- function(gc_df, bin_df) {
  a <- bin_df$bincount + 1
  bin_df$ratio <- a / mean(a)
  bin_df$lowratio <- lowess_gc(gc_df$gc.content, bin_df$ratio)
  
  return(bin_df)
}

cbs_segment_varbin <- function(bin_mat, alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5){
  assertthat::assert_that(is.numeric(bin_mat$chrom.numeric))
  assertthat::assert_that(all(!is.na(bin_mat$chrom.numeric)))

  set.seed(25)
  
  CNA.object <- DNAcopy::CNA(
      log(bin_mat$lowratio, base=2), 
      bin_mat$chrom.numeric,
      bin_mat$chrompos, data.type="logratio")
  smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- DNAcopy::segment(
      smoothed.CNA.object, alpha=alpha,
      nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2)
  #A table contains segmented regions and the corresponding segment means (log)
  thisShort <- segment.smoothed.CNA.object[[2]]
  
  #A long table:each row corresponds to the segment mean for each bin
  m <- matrix(data=0, nrow=nrow(bin_mat), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
    prevEnd = thisEnd
  }
  cbs.long.nobad <- m[, 1]
  
  #####  NEW STUFF  also check min.width=2 above
  
  workShort <- thisShort
  workShort$segnum <- 0
  workShort$seg.start <- 0
  workShort$seg.end <- 0
  
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    workShort$seg.start[i] <- thisStart
    workShort$seg.end[i] <- thisEnd
    workShort$segnum[i] <- i
    prevEnd = thisEnd
  }
  
  ### deal with short segments (appending, combining)
  discardSegments <- TRUE
  while (discardSegments) {
    orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
    if (orderShort[1, "num.mark"] < min.width) {
      workShort <- remove.segment(workShort, orderShort[1, "segnum"], bin_mat, undo.SD)
    } else {
      discardSegments <- FALSE
    }
  }
  
  workShort <- sdundo.all(workShort, bin_mat, undo.SD)
  thisShort <- workShort
  
  #####  END NEW STUFF
  
  m <- matrix(data=0, nrow=nrow(bin_mat), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
    prevEnd = thisEnd
  }
  
  bin_mat$seg.mean.LOWESS <- m[, 1]
  
  thisGrid <- seq(1.5, 5.5, by=0.05)
  thisOuter <- bin_mat$seg.mean.LOWESS %o% thisGrid
  thisOuterRound <- round(thisOuter)
  thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
  thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
  thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
  
  bin_mat$seg.quantal <- bin_mat$seg.mean.LOWESS * thisMultiplier
  bin_mat$ratio.quantal <- bin_mat$lowratio * thisMultiplier
  
  return(bin_mat)
}

