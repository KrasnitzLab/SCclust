# 
# Author: lubo
# Created: Mar 21, 2018
###############################################################################
library("assertthat")

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

cbs_segment_varbin <- function(
    bin_mat, alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5,
    mult.min=1.5, mult.max=5.5){

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
    for(i in 1:nrow(orderShort)) {
      if(sum(workShort[workShort[, "chrom"] == orderShort[i, "chrom"], "num.mark"]) >= min.width) {
        if(orderShort[i, "num.mark"] < min.width) {
          workShort <- remove_segment(workShort, orderShort[i, "segnum"], bin_mat, undo.SD)
          break
        } 
      }
      discardSegments <- FALSE
    }
  }
  
  workShort <- sdundo_all(workShort, bin_mat, undo.SD)
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
  
  thisGrid <- seq(mult.min, mult.max, by=0.05)
  thisOuter <- bin_mat$seg.mean.LOWESS %o% thisGrid
  thisOuterRound <- round(thisOuter)
  thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
  thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
  thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
  
  bin_mat$seg.quantal <- bin_mat$seg.mean.LOWESS * thisMultiplier
  bin_mat$ratio.quantal <- bin_mat$lowratio * thisMultiplier
  
  return(bin_mat)
}


sdundo_all <- function (sdShort, ratioData, sd.undo) {
  
  tempShort <- sdShort
  thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
  
  while ( TRUE ) {
    
    chrom <- tempShort$chrom
    chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])
    breakpoints <- which(chrom == chrom.shift)
    #breakpoints intra(inside) chrom
    
    if (length(breakpoints) < 1) {
      break
    }
    
    breakpoints.shift <- breakpoints + 1
    undo.breakpoints <- breakpoints[
        which(abs(tempShort$seg.mean[breakpoints] - 
                    tempShort$seg.mean[breakpoints.shift]) < thisSd * sd.undo)]
    
    if (length(undo.breakpoints) < 1) {
      break
    }
    
    undo.breakpoints.shift <- undo.breakpoints + 1
    undo.df <- tempShort[undo.breakpoints, ]
    undo.df$seg.mean.diff <- abs(
        tempShort$seg.mean[undo.breakpoints] - 
            tempShort$seg.mean[undo.breakpoints.shift])
    min.index <- which.min(undo.df$seg.mean.diff)
    leftIndex <- undo.df$segnum[min.index]
    rightIndex <- leftIndex + 1
    tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
    tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
    tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + 
        tempShort[rightIndex, "num.mark"]
    tempShort[leftIndex, "seg.mean"] <- mean(
        log(ratioData$lowratio[
                tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], 
            base=2))
    tempShort <- tempShort[-rightIndex, ]
    tempShort$segnum <- seq(1:nrow(tempShort))
  }
  
  return(tempShort)
}


remove_segment <- function( rsShort, rsSegnum, ratioData, sd.undo ) {
  
  ############################################
  ## deciding the appending location (left or right)
  ############################################
  
  appendLeft <- TRUE
  checkSdundo <- FALSE
  
  if (rsSegnum == 1) {
    appendLeft <- FALSE
  } else {
    if (rsSegnum == nrow(rsShort)) {
      appendLeft <- TRUE
    } else {
      rightIndex <- rsSegnum + 1
      leftIndex <- rsSegnum - 1
      
      if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
        appendLeft <- TRUE
      } else {
        if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
          appendLeft <- FALSE
        } else {
          if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) < 
              abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
            appendLeft <- TRUE
            checkSdundo <- TRUE
          } else {
            appendLeft <- FALSE
            checkSdundo <- TRUE
          }
        }
      }
    }
  }
  
  appendIndex <- 99999999
  if (appendLeft) {
    appendIndex <- rsSegnum - 1
  } else {
    appendIndex <- rsSegnum + 1
  }
  
  ## in the short table(each row is one segment), append the short segment 
  ## to its left/right segment.
  
  ############################################
  ## make the new short table after appending the short segment.
  ############################################
  
  tempShort <- rsShort
  newLocStart <- -1
  newLocEnd <- -1
  if (appendLeft) {
    tempShort[appendIndex, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
    tempShort[appendIndex, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
  } else {
    tempShort[appendIndex, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
    tempShort[appendIndex, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
  }
  
  tempShort[appendIndex, "num.mark"] <- tempShort[appendIndex, "num.mark"] + 
      tempShort[rsSegnum, "num.mark"]
  ## updating seg mean
  tempShort[appendIndex, "seg.mean"] <- mean(
      log(ratioData$lowratio[
              tempShort[appendIndex, "seg.start"]:tempShort[appendIndex, "seg.end"]], 
          base=2))
  
  tempShort <- tempShort[-rsSegnum, ]
  tempShort$segnum <- seq(1:nrow(tempShort))
  
  if (checkSdundo) {
    thisSd <- -1
    if (appendLeft) {
      leftIndex <- appendIndex
      rightIndex <- appendIndex + 1
    } else {
      leftIndex <- appendIndex - 2
      rightIndex <- appendIndex - 1
    }
    # sd of log2 (reason to use sqrt(2) here)GC normalized bin count data 
    # (the data before segmentation)
    thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
    
    ############################################
    ## after appending, check the new splitting.
    ############################################
    
    ## check whether the new splitting a real change point
    
    if (abs(tempShort$seg.mean[leftIndex] - tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {
      
      ##  remove changepoint (combine)
      tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
      tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
      tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + 
          tempShort[rightIndex, "num.mark"]
      ## updating seg mean
      tempShort[leftIndex, "seg.mean"] <- mean(
          log(ratioData$lowratio[
                  tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], 
              base=2))
      tempShort <- tempShort[-rightIndex, ]
      tempShort$segnum <- seq(1:nrow(tempShort))
    }
  }
  
  return(tempShort)
}

#' Generate the segmented profile for each cell.
#'
#' Generate the segmented profile for each cell in the input directory using CBS. 
#'
#' @param varbin_files list of bin count files for all cells produced by 
#'        'varbin' step of 'sgains' package.
#' @param gc_df binning scheme used for the analysis.
#' @param badbins list of bins that should be excluded from the analysis.
#' @return The list containing seg quantal and ratio quantal matrix for all cells.
#' @export
segment_varbin_files <- function(varbin_files, gc_df, badbins=NULL) {
  
  processed <- vector("list", nrow(varbin_files))
  
  ncells <- nrow(varbin_files)
  cells <- varbin_files$cells
  
  for(index in seq(1, ncells)) {
      varbin_file <- varbin_files$paths[index]
      cell <- varbin_files$cells[index]
      name <- varbin_files$names[index]
            
      bin_df <- load_table(varbin_file)
      bin_df$chrom.numeric <- chrom_numeric(bin_df$chrom)
      if(!is.null(badbins)) {
        bin_df <- bin_df[-badbins,]
      }
      assertthat::assert_that(nrow(bin_df) == nrow(gc_df))

      bin_df <- cbs_segment_ratio(gc_df, bin_df)
      bin_df <- cbs_segment_varbin(bin_df)
    
      processed[[index]] <- bin_df
  }
  
  uber_seg <- data.frame(
          chrom=gc_df$chrom.numeric, 
          chrompos=gc_df$bin.start, 
          abspos=gc_df$bin.start.abspos)
  
  uber_ratio <- data.frame(
          chrom=gc_df$chrom.numeric, 
          chrompos=gc_df$bin.start, 
          abspos=gc_df$bin.start.abspos)
  
  for(index in seq(1, ncells)) {
      bin_df <- processed[[index]]
      assertthat::assert_that(!is.null(bin_df$chrom))
      assertthat::assert_that(!is.null(bin_df$chrom.numeric))
      assertthat::assert_that(!is.null(bin_df$abspos))
      assertthat::assert_that(!is.null(bin_df$ratio.quantal))
      assertthat::assert_that(!is.null(bin_df$seg.quantal))

      cell <- cells[index]
      uber_seg[[cell]] <- bin_df$seg.quantal
      uber_ratio[[cell]] <- bin_df$ratio.quantal
  }
  
  return(list(seg=uber_seg,ratio=uber_ratio))
}