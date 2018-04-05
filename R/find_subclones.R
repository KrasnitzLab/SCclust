#'Identify subclones in hierarchical tree.
#'
#'Based on hierarchical clustering, identify the hard/soft clones.
#'@param hc The hclust objects with clones identified.
#'@param pinmat The pinmat.
#'@param pins The pins.
#'@param minnodesize An integer. Default: 6. The minimum node size for a subclone.
#'@param nsim The number of permutation simulations for subclone identification. Default: 500.
#'@param lmmax Numeric value. Default: 0.001. The threshold parameter for the linear fit to identify subclones.
#'@param hcmethod Default: average
#'@param base_share An integer. Default: 3. A balance parameter for controlling minimal number of shared features in a subclone node.
#'@param fdr_thresh FDR criterion for subclone nodes. Default: -2.
#'@param sharemin A feature is considered shared if present in share_min fraction of leaves in a node.Default: 0.90.
#'@param bymax Logical. If TRUE (Default), use maximal of mean FDR for the node to find subclones.
#'@param climbfromsize An integer. Default: 2.
#'@param climbtoshare An integer. Default: 3.
#'@return A list of hclust objects for clones.
#'@export


find_subclones <- function(
    hc, pinmat, pins, 
    nmin = -6, nsim = 500, lmmax = 0.001, 
    hcmethod = "average", baseshare = 3, 
    fdrthresh = -2, sharemin = 0.85, 
    bymax = TRUE, climbfromsize = 2, climbtoshare = 3, clonetype = 'soft'){

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

      n_share <- baseshare + sum(rowSums(clonemat) > (sharemin * ncol(clonemat)))
      cellnames <- colnames(clonemat)
      print(paste("subclone number of cells: ", length(cellnames)))
      res <- sim_fisher_wrapper(clonemat, clonepins, nsim=nsim)
      true_pv <- res$true
      assertthat::assert_that(all(!is.null(true_pv)))
      
      sim_pv <- res$sim
      sim_pv <- sim_pv[!is.na(sim_pv)]
      assertthat::assert_that(all(!is.null(sim_pv)))

      mfdr <- fisher_fdr(true_pv, sim_pv, cellnames, lmmax = lmmax)

      mdist <- fisher_dist(true_pv, cellnames)

      subhc <- hclust_tree(clonemat, mfdr, mdist, hcmethod = hcmethod)
      subhc_clone <- find_clones(
          subhc, fdrthresh = fdrthresh, sharemin = sharemin, 
          nshare = n_share, bymax = bymax,
          climbfromsize = climbfromsize, climbtoshare = climbtoshare)
      subhc_clones[[i]] <- subhc_clone
    }

  }

  return(construct_clonetable(hc, subhc_clones, clonetype))
}

construct_clonetable <- function(
    hc, subhc_clones, clonetype, subclone_toobig=0.8){
  assertthat::assert_that(!is.null(hc$labels))
  assertthat::assert_that(!is.null(hc$softclones))
  assertthat::assert_that(!is.null(hc$leaflist))
  
  clonetable <- data.frame(hc$labels, rep(0, length(hc$labels)),
      rep(0, length(hc$labels)), stringsAsFactors = F)
  colnames(clonetable) <- c("ID", "clone", "subclone")
  
  ##clones
  for (nodes in unique(hc$softclones[clonetype,])){
    clonetable[hc$leaflist[[nodes]], "clone"] <- nodes
  }
  
  ##subclones
  for (subhc in subhc_clones){
    assertthat::assert_that(!is.null(subhc$labels))
    assertthat::assert_that(!is.null(subhc$softclones))
    assertthat::assert_that(!is.null(subhc$leaflist))
    assertthat::assert_that(!is.null(subhc$nodesize))
    
    clunique <- unique(subhc$softclones[clonetype,])
    print(clunique)

    if (length(clunique) > 1){
      for(nodes in clunique){
        clonetable[
            match(subhc$labels[[nodes]], clonetable[,1]), 
            "subclone"] <- nodes
      }
    } else if(length(clunique) == 1){
      if(subhc$nodesize[clunique] < (subclone_toobig * max(subhc$nodesize))){
        clonetable[
            match(subhc$labels[[clunique]], clonetable[,1]), 
            "subclone"] <- clunique}
    }
  }
  
  return(clonetable)
}


#find_subclones <- function(
#    hc, pinmat, pins, min_node_size = -6, sim_round = 500, lm_max = 0.001, 
#    hc_method = "average", base_share = 3, fdr_thresh = -2, share_min = 0.90, 
#    bymax = TRUE, climb_from_size = 2, climb_to_share = 3, clonedef = 'soft',
#    clonetype="soft", subcloneTooBig=0.8) {
#
#  clonetable<-data.frame(hc$labels,rep(0,length(hc$labels)),
#      rep(0,length(hc$labels)),stringsAsFactors=F)
#  
#  dimnames(clonetable)[[2]]<-c("ID","clone","subclone")
#  for(nodes in unique(hc$softclones[clonetype,])) {
#    clonetable[hc$leaflist[[nodes]],"clone"]<-nodes
#  }
#  
#  sub_hc_clones <- find_subclones_calc(
#      hc, pinmat, pins, min_node_size, sim_round, lm_max, 
#      hc_method, base_share, fdr_thresh, share_min, 
#      bymax, climb_from_size, climb_to_share, clonedef)
#  
#  for(hc_clone in sub_hc_clones) {
#    print(head(hc_clone))
#
#    clunique<-unique(hc_clone$softclones[clonetype,])
#    
#    if(length(clunique)>1) {
#      for(nodes in clunique){
#        clonetable[match(hc_clone$labellist[[nodes]],clonetable[,1]),"subclone"]<-nodes
#      }
#    } else if(length(clunique)==1) {
#      if(hc_clone$nodesize[clunique]<(subcloneTooBig*max(hc_clone$nodesize)))
#        clonetable[match(hc_clone$labellist[[clunique]],clonetable[,1]),"subclone"]<-
#            clunique
#    }
#  }
#  return(clonetable)
#}



