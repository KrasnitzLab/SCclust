#'Compute FDRs for Fisher's test p-values.
#'
#'Linear fit to the tail of empirical null distribution of Fisher p-values;
#'FDR computation: compare true to simulated CDF(empirical null).
#'@param true_pv The Fisher's test p-values for the observation.
#'@param sim_pv The Fisher's test p-values for the permutations.
#'@param cell_names A character vector. The names of cells.
#'@param lm_max Numeric value. Default: 0.001. The threshold parameter 
#'          for the linear fit.
#'@return A list containing the matrix of the FDR values (mat_fdr) and 
#'          distance matrix based on Fisher's test p-values (mat_dist).
#'@export


fisher_fdr <- function(
    true_pv, sim_pv, cell_names, lm_max = 0.001){

  assertthat::assert_that(length(cell_names) == (1 + sqrt(1 + 8*length(true_pv)))/2)
      
      
  # Sort the true and the null data and get counts for each unique 
  # fisher pv value
  sim_sort <- sort(sim_pv)
  sim_unique <- unique(sim_sort)
  sim_count <- tapply(
          match(sim_sort, sim_unique), 
          match(sim_sort, sim_unique), length)
  true_sort <- sort(true_pv)
  true_unique <- unique(true_sort)
  true_count <- tapply(
          match(true_sort, true_unique), 
          match(true_sort, true_unique), length)
  
  
  # linear fit to the tail of empirical null distribution of Fisher p-values
  ##############################################################################
  
  # lm_max -- A parameter in a linear fit
  # Observed p-values are often far lower than any p-value sampled from the null; 
  # to determine FDR in such cases get a power-law fit to the low-p tail of the 
  # null CDF and use it to extrapolate to very low p-values. Use the actual 
  # null CDF to estimate FDR for higher p-values.
  lowp_index <- (cumsum(sim_count)/sum(sim_count)) < lm_max
  lmfit <- lm(
          log(cumsum(
                sim_count[lowp_index])/sum(sim_count))~
          log(sim_unique[lowp_index]))
  
  ## FDR computation: compare true to simulated CDF(empirical null)
  ##############################################################################
  dat_TrueSim <- cbind(
      c(log(true_unique), log(sim_unique)),
      c(log(cumsum(true_count)/sum(true_count)), 
          log(cumsum(sim_count)/sum(sim_count))),
      c(rep(0, length(true_unique)), rep(1, length(sim_unique))))
  dat_TrueSim <- dat_TrueSim[order(dat_TrueSim[,1]),]
  
  ## (1) estimate FDR for higher p-values
  # (linear interpolation)

  # num of sim with less pv than the specific true value
  num_sim <- cumsum(dat_TrueSim[,3])[dat_TrueSim[,3] == 0]  
  x1pos <- match(num_sim, cumsum(dat_TrueSim[,3]))[num_sim > 0]
  x2pos <- match(num_sim + 1, cumsum(dat_TrueSim[,3]))[num_sim > 0]
  logpv1 <- dat_TrueSim[x1pos,1]
  logpv2 <- dat_TrueSim[x2pos,1]
  logcdf1 <- dat_TrueSim[x1pos,2]
  logcdf2 <- dat_TrueSim[x2pos,2]
  
  logfdr <- rep(0, length(true_unique))
  logfdr[num_sim > 0] <- 
      (logcdf2 - logcdf1)*log(true_unique)[num_sim > 0]/(logpv2 - logpv1) +
      (logcdf1*logpv2 - logcdf2*logpv1)/(logpv2 - logpv1) - 
      log(cumsum(true_count)/sum(true_count))[num_sim > 0]
  logfdr[true_unique > max(sim_unique)] <- 0
  
  # (2) estimate FDR for low p-values
  # (a power-law fit to the low-p tail of the null CDF, use it to extrapolate 
  # to very low p-values)
  lmu <- max(sim_unique[lowp_index])
  
  logfdr_lmfit <-
          lmfit$coefficients[2]*log(true_unique) + 
          lmfit$coefficients[1] - log(cumsum(true_count)/sum(true_count))
  
  if(is.finite(lmu) & min(true_unique) < min(sim_unique)){
    interp_index <- true_unique < lmu & true_unique > min(sim_unique)
    logfdr[interp_index]<-
        (logfdr_lmfit[interp_index] * (log(lmu) - log(true_unique[interp_index])) -
          logfdr[interp_index]*(log(min(sim_unique)) -
            log(true_unique[interp_index])))/(log(lmu) - log(min(sim_unique)))

    nointerp_index <- true_unique < min(sim_unique)
    logfdr[nointerp_index] <- logfdr_lmfit[nointerp_index]
  }
  
  logfdr <- cummax(logfdr)
  logfdr[logfdr > 0] <- 0
  
  logfdr_all <- logfdr[match(true_pv, true_unique)]
#  size <- length(cell_names)
#  mat_fdr <- matrix(ncol = size, nrow = size, data = 0)
#  mat_fdr[upper.tri(mat_fdr)] <- logfdr_all
#  mat_fdr <- pmin(mat_fdr,t(mat_fdr))
#    
#  colnames(mat_fdr) = rownames(mat_fdr) <- cell_names
  mat_fdr <- cell2cell_matrix(logfdr_all/log(10), cell_names)  # ??log(10)
  
#  mat_dist <- matrix(ncol = size, nrow = size, data = 0)
#  mat_dist[upper.tri(mat_dist)] <- log10(true_pv)
#  mat_dist <- pmin(mat_dist, t(mat_dist))
#  colnames(mat_dist) = rownames(mat_dist) <- cell_names

  mat_dist <- cell2cell_matrix(log10(true_pv), cell_names)
  
  return(list(mat_fdr = mat_fdr, mat_dist = mat_dist))
}


cell2cell_matrix <- function(data, cell_names) {
  size <- length(cell_names)
  assertthat::assert_that(length(cell_names) == (1 + sqrt(1 + 8*length(data)))/2)
  
  res <- matrix(ncol = size, nrow = size, data = 0)
  res[upper.tri(res)] <- data
  res <- pmin(res,t(res))
  
  colnames(res) = rownames(res) <- cell_names
  return(res)
}