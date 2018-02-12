

calc_centroareas <- function(cyto) {
  centromere<-c("p11","q11")
  cyto[,1]<-numeric_chrom(cyto[,1])
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

calc_badbins <- function(gc, centroareas) {
  # print(head(gc))

  gc$chrom <- numeric_chrom(gc$bin.chrom)


  for(i in centroareas$chrom) {
    from <- centroareas[i,]$from
    to <- centroareas[i, ]$to
    # print(from)
    # print(to)

    df <- gc[gc$chrom == i, ]
    bad_df <- df[df$bin.start >= from && df$bin.end <= to, ]
    # print(bad_df)

  }
}
