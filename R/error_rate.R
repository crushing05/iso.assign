#' Estimate the assignment error rate for known-origin tissue samples
#'
#' Function to estimate the proportion of individuals that are incorrectly assigned to their actual breeding origin
#' @param summ Assignment results from either the iso_assign() or abun_assign() functions
#' @param origin_cell Vector containing the true breeding cell for each individual
#' @param iso Should area be estimated for the isotope-only assignment (default) or with abundance used a prior
#' @keywords assignment; stable isotopes; error rate
#' @return Error rate (i.e., the proportion of individuals that are miss classified)
#' @examples
#' error_rate(summ = woth_iso, origin.cell = woth_origin_cell)
#'

error_rate <- function(summ, origin_cell, iso = TRUE){

  if(iso == TRUE){
    # Convert vector of likely/unlikely classifications to cell x indv matrix
    origin <- matrix(summ$iso_origin, ncol = length(unique(summ$indv)), byrow = FALSE)
    # Extract classification for the true origin of each individual
    correct <- origin[cbind(origin_cell, seq_along(origin_cell))]
    # Estimate error rate
    err_rate <- 1 - mean(correct)
  }else{
    origin <- matrix(summ$wght_origin, ncol = length(unique(summ$indv)), byrow = FALSE)
    correct <- origin[cbind(origin_cell, seq_along(origin_cell))]
    err_rate <- 1 - mean(correct)
  }

  return(err_rate)
}
