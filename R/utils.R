library("futile.logger")

filter_evil_short <- function(short_df, eviltwins=NULL) {
  good_cells <- unique(short_df[,"profid"])
  if(!is.null(eviltwins)) {
    good_cells <- setdiff(good_cells, eviltwins)
    short_df <- short_df[short_df[, "profid"] %in% good_cells, ]
  }
  return(short_df)
}


filter_good_short <- function(short_df, good_cells) {
  short_df <- short_df[short_df[, "profid"] %in% good_cells, ]
  return(short_df)
}


filter_good_columns <- function(df, good_cells, skip=0) {
  cols <- colnames(df) %in% good_cells
  cols[1:skip] <- TRUE
  return(df[, cols])
}


filter_evil_columns <- function(df, evil_cells, skip=0) {
  cols <- ! colnames(df) %in% evil_cells
  cols[1:skip] <- TRUE
  return(df[, cols])
}

#' Calculates centromere regions (areas). 
#' 
#' @export
calc_centroareas <- function(cyto) {
  centromere<-c("p11","q11")
  cyto[,1]<-chrom_numeric(cyto[,1])
  cyto<-cyto[order(cyto[,1]),]

  centroleft<-cyto[grep(centromere[1],cyto[,4]),]
  centroright<-cyto[grep(centromere[2],cyto[,4]),]
  centroleft<-centroleft[match(unique(centroleft[,1]),centroleft[,1]),]
  centroright<-centroright[nrow(centroright):1,]
  centroright<-centroright[match(unique(centroright[,1]),centroright[,1]),]
  centroright<-centroright[nrow(centroright):1,]
  centroareas<-cbind(centroleft[,c(1,2)],centroright[,3])
  dimnames(centroareas)[[2]]<-c("chrom","from","to")

  return(centroareas)
}

#' Converts regions to list of bins from binning scheme. 
#' 
#' @export
calc_regions2bins <- function(gc_df, regions) {

  assertthat::assert_that(!is.null(gc_df$chrom.numeric))

  bins <- list()
  for(index in seq(nrow(regions))) {
    region <- regions[index, ]

    from <- region$from
    to <- region$to
    chrom <- region$chrom

    df <- gc_df[gc_df$chrom.numeric == chrom,]
   
    df <- df[
        ((df$bin.start >= from) & (df$bin.start <= to)) | 
            ((df$bin.end >= from) & (df$bin.end <= to)),]

    flog.debug("chrom: %s; from: %s; to: %s; bins filtered: %s", 
        chrom, from, to, nrow(df))
    bins <- append(bins, rownames(df))
  }
  return(as.numeric(unlist(bins, recursive=T)))
}


#' Converts list of bins from binning scheme to regions. 
#' 
#' @export
calc_bins2regions <- function(gc_df, bins) {
  regions <- NULL
  region <- NULL

  for(bin in bins) {
    d <- gc_df[bin,]
    chrom <- d$chrom.numeric
    from <- d$bin.start
    to <- d$bin.end
    
    if(is.null(region)) {
      region <- data.frame(chrom, from, to)
    } else if(region$to == from & region$chrom == chrom) {
      # extend
      flog.debug("extending region: %s", region)
      region$to <- to
    } else {
      if(is.null(regions)) {
        regions <- region
      } else {
        regions <- rbind(regions, region)
      }
      region <- data.frame(chrom, from, to)
    }
  }
  regions <- rbind(regions, region)
  return(regions)
}

#' Converts chrom name to numeric and adds `chrom.numeric` column to the 
#' dataframe. 
#' 
#' @export
chrom_numeric <- function(chrom) {
  if(is.numeric(chrom)) {
    chrom.numeric <- chrom
    return(chrom.numeric)
  } else {
    chrom.numeric <- substring(chrom, 4)
    chrom.numeric[which(chrom == "chrX")] <- "23"
    chrom.numeric[which(chrom == "chrY")] <- "24"
    
    chrom.numeric <- as.numeric(chrom.numeric)
    return(chrom.numeric)
  }
}


chrom_numeric_mouse <- function(chrom) {
  if(is.numeric(chrom)) {
    chrom.numeric <- chrom
    return(chrom.numeric)
  } else {
    chrom.numeric <- substring(chrom, 4)
    chrom.numeric[which(chrom == "chrX")] <- "20"
    chrom.numeric[which(chrom == "chrY")] <- "21"
    
    chrom.numeric <- as.numeric(chrom.numeric)
    return(chrom.numeric)
  }
}


tree_clustersize <- function(indextable) {
  clustersize<-rep(NA,nrow(indextable))

  csleft<-rep(NA,nrow(indextable))
  csleft[indextable[,"index1"]<0]<-1

  csright<-rep(NA,nrow(indextable))
  csright[indextable[,"index2"]<0]<-1

  while(is.na(sum(clustersize))){
    clustersize<-csleft+csright
    csleft[indextable[,"index1"]>0]<-
        clustersize[indextable[indextable[,"index1"]>0,"index1"]]
    csright[indextable[,"index2"]>0]<-
        clustersize[indextable[indextable[,"index2"]>0,"index2"]]
  }  

  return(clustersize)
}


#' Builds HC tree representation based on the distance matrix
#' computed by \code{fisher_dist}
#' 
#' @export
tree_py <- function(mdist, method, metric='euclidean'){
  hc<-hclust(as.dist(mdist), method)
  
  res <- cbind(hc$merge, hc$height)
  colnames(res) <- c("index1", "index2", "height")
  clustersize <- tree_clustersize(res)
  
  d <- res[, 1:2]
  d[res[,1:2]<0]<- -d[res[,1:2]<0]-1
  d[res[,1:2]>0]<- d[res[,1:2]>0]+nrow(res)
  res[,1:2]<-d

  res <- cbind(res, clustersize)
  colnames(res)<-c("index1","index2","height", "clustersize")
  return(res)
}

hc_climb<-function(hc, minsize, minshare){
  hard2soft<-matrix(nrow=2,ncol=0,dimnames=list(c("hard","soft"),NULL))
  for(cnode in hc$clonenodes) {
    if(hc$nodesize[cnode]>minsize){
      nodenow<-cnode
      ancestor<-cnode
      while(nodenow<nrow(hc$merge)){
        nodenow<-row(hc$merge)[hc$merge==nodenow]
        if(hc$count_pins_share[nodenow]>=minshare)ancestor<-nodenow
      }
      
      hard2soft<-cbind(hard2soft,c(cnode,ancestor))
    }
  }
  return(hard2soft)
}

