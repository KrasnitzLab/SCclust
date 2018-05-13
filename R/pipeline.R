# 
# Author: lubo
# Created: Mar 28, 2018
###############################################################################
library("futile.logger")


#'Integration to 'sGAINS' tool. 
#' 
#' This function is called by sGAINS tools to perform final step in preparation
#' of results.
#'
#' @param scgv_dir directory where the results of the analysis should be stored
#' @param case_name name of the case to be used for storing results of the analysis
#' @param varbin_dir directory where are located results from 'varbin' step of sGAINS
#' @param varbin_suffix common suffix for files produced from 'varbin' step of sGAINS
#' @param bins_boundaries file name for binning scheme to use in the analysis
#' @param cytoband file name where is located the description of cyto bands for
#'        the version of genome we are using
#' @param badbins a file name were the defenition of bad bins for the specified
#'        binning schema could be found
#' @param nsim number of simulations to run for calculating simulated FDR distribution
#' @param sharemin a feature is considered shared if present in sharemin fraction 
#'        of leaves in a node
#' @export

sgains_pipeline <- function(
    scgv_dir, case_name,
    varbin_dir, varbin_suffix,
    bins_boundaries,
    cytoband,
    badbins=NULL,
    nsim=150,
    sharemin=0.85) {

  if(!file.exists(scgv_dir)) {
    dir.create(scgv_dir)
  }
  filenames <- case_filenames(scgv_dir, case_name)
  
  flog.debug("loading bin boundarires...")
  gc_df <-load_table(bins_boundaries)
  gc_df$chrom.numeric <- chrom_numeric(gc_df$bin.chrom)
  
  if(is.null(badbins)) {
    flog.debug("loading cytoband data...")
    cyto_data <- read.table(cytoband, header=F,as.is=T)
    flog.debug("calculating centromere regions...")
    dropareas <- calc_centroareas(cyto_data)
    flog.debug("calculating dropbins from centromere regions...")
    badbins <- calc_regions2bins(gc_df, dropareas)
  }
  
  gc_df <- gc_df[-badbins, ]
  flog.debug("centrobins filtered; gc rows: %s", nrow(gc_df))
  
  flog.debug("looking for varbin files...")
  varbin_files <- varbin_input_files(varbin_dir, varbin_suffix)
  cells <- varbin_files$cells
  
  flog.debug("varbin files found: ", varbin_files)
  print(head(varbin_files))
  
  res <- segment_varbin_files(varbin_files, gc_df, badbins)
  uber_seg <- res$seg
  uber_ratio <- res$ratio
  
  cells <- uber_cells(uber_seg, skip=3)$cells

  save_table(filenames$cells, data.frame(cell=cells))
  save_table(filenames$seg, uber_seg)
  save_table(filenames$ratio, uber_ratio)
  
  pins <- calc_pinmat(gc_df, uber_seg, dropareas=dropareas)
  pinmat_df <- pins$pinmat
  pins_df <- pins$pins
  save_table(filenames$featuremat, pinmat_df)
  save_table(filenames$features, pins_df)
  
  fisher <- sim_fisher_wrapper(
      pinmat_df, pins_df, njobs=30, nsim=nsim, nsweep=10)
  true_pv <- fisher$true
  sim_pv <- fisher$sim
  
  mfdr <- fisher_fdr(true_pv, sim_pv, cells)
  mdist <- fisher_dist(true_pv, cells)
  
  hc <- hclust_tree(pinmat_df, mfdr, mdist)
  
  tree_df <- tree_py(mdist, method='average')

  hc <- find_clones(hc, sharemin=sharemin)
  subclones <- find_subclones(hc, pinmat_df, pins_df, nsim=nsim, sharemin=sharemin)
  
  save_table(filenames$tree, tree_df)
  save_table(filenames$clone, subclones)
  
}