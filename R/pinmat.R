

augment_gc <- function(gc_df, df) {
  assertthat::assert_that(all(gc_df$chrom.numeric == df$chrom))
  assertthat::assert_that(all(gc_df$bin.start == df$chrompos))

  augment_df <- cbind(gc_df[,c("chrom.numeric", "bin.start", "bin.end")], df[,"abspos"])
  colnames(augment_df) <- c("chrom", "chromstart", "chromend", "absstart")

  augment_df <- cbind(
    augment_df,
    augment_df[,"absstart"] + (augment_df[,"chromend"] - augment_df[,"chromstart"])
  )
  colnames(augment_df) <- c("chrom", "chromstart", "chromend", "absstart", "absend")
  return(augment_df)
}

calc_segments_short <- function(gc_df, segment_df, homoloss=0.01) {
  assertthat::assert_that(is.numeric(gc_df$chrom))
  assertthat::assert_that(is.numeric(segment_df$chrom))
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

  tshort <- cbind(tshort, tshort[,"segvals"] - 
          ploidies_df[tshort[,"profid"], "ploidychromod"])
  colnames(tshort)[ncol(tshort)] <- "cvals"

  if(homoloss > 0) {
    tshort <- tshort[tshort[,"profid"]%in%dimnames(ploidies_df)[[1]][
            ploidies_df[,"homoloss"] <= homoloss],]
  }
  return(tshort)
}

calc_ploidies <- function(gc_df, segment_df, a=NULL) {
  if(is.null(a)) {
    a <- round(as.matrix(segment_df[,-(1:3)]))
  }
  assertthat::assert_that(nrow(a) == nrow(segment_df))
  assertthat::assert_that(ncol(a) + 3 == ncol(segment_df))

  getmode<-function(x) as.numeric(
        names(which.max(tapply(X=x,INDEX=as.factor(x), FUN=length))))
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

calc_censored_index <- function(short_df, dropareas) {
  assertthat::assert_that(!is.null(dropareas))
  censored <-
    ((short_df[, "chromstart"] >= dropareas[short_df[, "chrom"], "from"]) &
       (short_df[, "chromend"] <= dropareas[short_df[, "chrom"], "to"]))
  censoredtoo <-
    (which(censored) + 1)[(which(censored) + 1) <= nrow(short_df)]

  censored[censoredtoo] <-
    ((short_df[censoredtoo, "chrom"] == short_df[censoredtoo - 1, "chrom"]) &
       (short_df[censoredtoo, "profid"] == short_df[censoredtoo - 1, "profid"])) | 
     censored[censoredtoo]
  return(censored)
}

filter_dropareas_short <- function(short_df, dropareas = NULL) {
  if (!is.null(dropareas)) {
    censored <- calc_censored_index(short_df, dropareas)
    return(short_df[!censored,])
  }
  return(short_df)
}


calc_smear_breakpoints <- function(
    short_df, censored, smear=1, keepboundaries=F, chromrange=1:22) {

  dtshort<-cbind(
    short_df[,c("profid","chrom")],
    short_df[,"segstarts"]-smear,
    short_df[,"segstarts"]+smear,
    sign(short_df[,"cvals"]-c(0,short_df[-nrow(short_df),"cvals"]))
  )
  dimnames(dtshort)[[2]]<-c("profid","chrom","bpstart","bpend","bpsign")

  ustart<-short_df[match(unique(short_df[,"chrom"]),short_df[,"chrom"]),"segstarts"]
  uend<-c((ustart-1)[-1],short_df[nrow(short_df),"segends"])

  if(keepboundaries){
    dtshort[(dtshort[,"bpstart"]+smear)==ustart[dtshort[,"chrom"]],"bpsign"]<-
      2*dtshort[(dtshort[,"bpstart"]+smear)==ustart[dtshort[,"chrom"]],"bpsign"]
    dtshort<-dtshort[!censored,]
  } else {
    dtshort<-dtshort[((dtshort[,"bpstart"]+smear) > ustart[dtshort[,"chrom"]])&!censored,]
    dtshort[dtshort[,"bpstart"]<ustart[dtshort[,"chrom"]],"bpstart"]<-
      ustart[dtshort[dtshort[,"bpstart"]<ustart[dtshort[,"chrom"]],"chrom"]]
    dtshort[dtshort[,"bpend"]>uend[dtshort[,"chrom"]],"bpend"]<-
      uend[dtshort[dtshort[,"bpend"]>uend[dtshort[,"chrom"]],"chrom"]]
  }
  dtshort<-dtshort[dtshort[,"chrom"]%in%chromrange,]
  return(dtshort)
}

containment.indicator <- function(vstart, vend, wstart, wend){
  lw<-length(wstart)
  lv<-length(vstart)
  z<-cbind(c(vend,wend),c(1:lv,rep(0,lw)),c(rep(0,lv),1:lw))
  z<-z[order(z[,1]),]
  endbeforeend<-cummax(z[,2])[order(z[,3])][sort(z[,3])!=0]
  z<-cbind(c(wstart,vstart),c(rep((lv+1),lw),1:lv),c(1:lw,rep(0,lv)))
  z<-z[order(z[,1]),]
  startafterstart<-rev(cummin(rev(z[,2])))[order(z[,3])][sort(z[,3])!=0]
  return(cbind(startafterstart,endbeforeend))
}

calc_pinmat_short <- function(short_df, smear_df) {
  allsigns<-sort(unique(smear_df[,"bpsign"]))
  for(vsign in allsigns){
    a<-smear_df[smear_df[,"bpsign"]==vsign,]
    a<-a[order(a[,"bpend"]),]
    apins<-NULL
    while(nrow(a)>0){
      apins<-c(apins,a[1,"bpend"])
      a<-a[!(a[,"bpstart"]<=apins[length(apins)] & 
                a[,"bpend"]>=apins[length(apins)]),,drop=F]
    }
    a<-smear_df[smear_df[,"bpsign"]==vsign,]
    ci<-containment.indicator(
        apins,
        apins,
        a[order(a[,"bpend"]),"bpstart"],
        a[order(a[,"bpend"]),"bpend"])
    a<-cbind(a[order(a[,"bpend"]),],ci)

    dimnames(a)[[2]][(ncol(a)-1):ncol(a)]<-c("startpin","endpin")
    apinmat<-matrix(
        ncol=length(unique(short_df[,"profid"])),
        nrow=length(apins)+2,
        data=0,
        dimnames=list(NULL,unique(short_df[,"profid"])))

    for(id in unique(a[,"profid"])){
      apinmat[a[a[,"profid"]==id,"startpin"]+1,id]<-1
      apinmat[a[a[,"profid"]==id,"endpin"]+2,id] <- 
          apinmat[a[a[,"profid"]==id,"endpin"]+2,id]-1
    }
    apinmat<-apply(apinmat,2,cumsum)
    apinmat<-apinmat[-c(1,nrow(apinmat)),,drop=F]
    apins<-cbind(apins,rep(vsign,length(apins)))
    dimnames(apins)[[2]]<-c("bin","sign")
    if(vsign==min(allsigns)){
      pinmat<-apinmat
      pins<-apins
    }
    else{
      pinmat<-rbind(pinmat,apinmat)
      pins<-rbind(pins,apins)
    }
  }
  pins<-pins[(rowSums(pinmat)<ncol(pinmat)),,drop=F]
  pinmat<-pinmat[(rowSums(pinmat)<ncol(pinmat)),,drop=F]
  pinmat_df <- data.frame(pinmat)
  pins_df <- data.frame(pins)

  res <- list(pinmat_df, pins_df)
  names(res) <- c("pinmat", "pins")

  return(res)
}


calc_pinmat <- function(gc_df, segment_df, homoloss=0.0, dropareas=NULL) {

  augment_df <- augment_gc(gc_df, segment_df)
  short_df <- calc_segments_short(augment_df, segment_df, homoloss=homoloss)
  
  censored_index <- calc_censored_index(short_df, dropareas)
  smear_df <- calc_smear_breakpoints(short_df, censored_index)
  
  return(calc_pinmat_short(short_df, smear_df))
}
