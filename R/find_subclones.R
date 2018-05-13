#' Identify subclones in hierarchical tree.
#'
#' Based on hierarchical clustering, identify the hard/soft clones.
#' @param hc The hclust objects with clones identified.
#' @param pinmat The pinmat.
#' @param pins The pins.
#' @param minnodesize An integer. Default: 6. The minimum node size for a subclone.
#' @param nsim The number of permutation simulations for subclone identification. 
#'        Default: 500.
#' @param lmmax Numeric value. Default: 0.001. The threshold parameter for the 
#'        linear fit to identify subclones.
#' @param hcmethod Default: average
#' @param baseshare An integer. Default: 3. A balance parameter for controlling 
#'        minimal number of shared features in a subclone node.
#' @param fdr_thresh FDR criterion for subclone nodes. Default: -2.
#' @param sharemin A feature is considered shared if present in sharemin fraction 
#'        of leaves in a node.Default: 0.85.
#' @param bymax Logical. If TRUE (Default), use maximal of mean FDR for the node 
#'        to find subclones.
#' @param climbfromsize An integer. Default: 2.
#' @param climbtoshare An integer. Default: 3.
#' @return A list of hclust objects for clones.
#' @export


find_subclones <- function(
    hc, pinmat, pins, 
    nmin = -6, 
    nsim = 500, 
    lmmax = 0.001, 
    hcmethod = "average", 
    baseshare = 3, 
    fdrthresh = -2, 
    sharemin = 0.85, 
    bymax = T, 
    climbfromsize = 2, 
    climbtoshare = 3, 
    clonetype = 'soft'){

  assertthat::assert_that(!is.null(hc$softclones))
  assertthat::assert_that(!is.null(hc$leaflist))
  assertthat::assert_that(!is.null(hc$labellist))
  assertthat::assert_that(!is.null(hc$labels))
  
  
  bigclones<- unique(hc$softclones[clonetype,])[
      hc$nodesize[unique(hc$softclones[clonetype,])] >= nmin]
  
  subhc_clones <- list()

  if(length(bigclones) > 0){

    for (i in 1:length(bigclones)){
      clonemat <- pinmat[,dimnames(pinmat)[[2]] %in%
              hc$labellist[[bigclones[i]]], drop = F]
      clonepins <- pins[
          (rowSums(clonemat) > 0) & (rowSums(clonemat) < ncol(clonemat)),,drop = F]
      clonemat <- clonemat[
          (rowSums(clonemat) > 0)&(rowSums(clonemat) < ncol(clonemat)),,drop = F]

      nshare <- baseshare + sum(rowSums(clonemat) > (sharemin * ncol(clonemat)))
      cellnames <- colnames(clonemat)
      flog.debug("subclone number of cells: %s", length(cellnames))
      res <- sim_fisher_wrapper(clonemat, clonepins, nsim=nsim)
      if(is.null(res)) {
        next
      }

      true_pv <- res$true
      assertthat::assert_that(all(!is.null(true_pv)))
      
      sim_pv <- res$sim
      sim_pv <- sim_pv[!is.na(sim_pv)]
      assertthat::assert_that(all(!is.null(sim_pv)))

      mfdr <- fisher_fdr(true_pv, sim_pv, cellnames, lmmax = lmmax)

      mdist <- fisher_dist(true_pv, cellnames)

      subhc <- hclust_tree(pinmat, mfdr, mdist, hcmethod = hcmethod)
      flog.debug("fdrthresh=%s; sharemin=%s; nshare=%s; bymax=%s;", 
          fdrthresh, sharemin, nshare, bymax)
      flog.debug("climbfromsize=%s; climbtoshare=%s",
          climbfromsize, climbtoshare)
      flog.debug("subhc=%s", subhc)

      subhc <- find_clones(
          subhc, 
          fdrthresh = fdrthresh, 
          sharemin = sharemin, 
          nshare = nshare, 
          bymax = bymax,
          climbfromsize = climbfromsize, 
          climbtoshare = climbtoshare)
      subhc_clones[[i]] <- subhc
    }

  }

  return(construct_clonetracks(hc, subhc_clones, clonetype))
}

construct_clonetracks<- function(
    hc, subhc_clones, clonetype, subclone_toobig=0.8){
  assertthat::assert_that(!is.null(hc$labels))
  assertthat::assert_that(!is.null(hc$softclones))
  assertthat::assert_that(!is.null(hc$leaflist))
  assertthat::assert_that(!is.null(hc$labellist))
  
  clonetable <- data.frame(
      hc$labels, 
      rep(0, length(hc$labels)),
      rep(0, length(hc$labels)), stringsAsFactors = F)
  colnames(clonetable) <- c("ID", "clone", "subclone")
  
  ##clones
  for (nodes in unique(hc$softclones[clonetype,])){
    clonetable[hc$leaflist[[nodes]], "clone"] <- nodes
  }
  
  ##subclones
  for (subhc in subhc_clones){
    if(is.null(subhc$labels)) {
      next
    }
    assertthat::assert_that(!is.null(subhc$labels))
    assertthat::assert_that(!is.null(subhc$softclones))
    assertthat::assert_that(!is.null(subhc$leaflist))
    assertthat::assert_that(!is.null(subhc$nodesize))
    assertthat::assert_that(!is.null(subhc$labellist))
    
    clunique <- unique(subhc$softclones[clonetype,])
    flog.debug("clunique=%s", clunique)

    if (length(clunique) > 1){
      for(nodes in clunique){
        clonetable[
            match(subhc$labellist[[nodes]], clonetable[,1]), 
            "subclone"] <- nodes
      }
    } else if(length(clunique) == 1){
      if(subhc$nodesize[clunique] < (subclone_toobig * max(subhc$nodesize))){
        clonetable[
            match(subhc$labellist[[clunique]], clonetable[,1]), 
            "subclone"] <- clunique}
    }
  }
  
  return(clonetable)
}




