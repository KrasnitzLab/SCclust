

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

calc_centrobins <- function(gc, centroareas) {

  assertthat::assert_that(!is.null(gc$chrom.numeric))

  badbins <- list()
  for(i in centroareas$chrom) {
    from <- centroareas[i,]$from
    to <- centroareas[i, ]$to

    bad_df <- gc[gc$chrom.numeric == i,]
   
    bad_df <- bad_df[
        ((bad_df$bin.start >= from) & (bad_df$bin.start <= to)) | 
            ((bad_df$bin.end >= from) & (bad_df$bin.end <= to)),]

    flog.debug("chrom: %s; from: %s; to: %s; bins filtered: %s", 
        i, from, to, nrow(bad_df))
    badbins <- append(badbins, rownames(bad_df))
  }
  return(as.numeric(unlist(badbins, recursive=T)))
}


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
