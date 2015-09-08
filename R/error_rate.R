#' Estimate the assignment error rate for known-origin tissue samples
#'
#' Function to estimate the proportion of individuals that are incorrectly assigned to their actual breeding origin
#' @param origin Matrix with columns indicating likely/unlikely cells for each individual.
#' @param origin.cell  Vector indicating the true breeding cell (i.e. row in origin matrix) for each individual.
#' @keywords assignment; stable isotopes; error rate
#' @return Error rate (i.e., the proportion of individuals that are miss classified)
#' @examples
#' error_rate(origin = amre_assign$iso.origin, origin.cell = amre$origin)
#'

error_rate <- function(origin, origin.cell){
  correct <- integer()
  for(i in 1:length(origin.cell)){
    correct[i] <- origin[origin.cell[i],i]==1
  }
  err_rate <- 1- mean(correct)

  return(err_rate)
}
