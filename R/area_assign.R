#' Estimate the assignment area for known-origin tissue samples
#'
#' Function to estimate the mean number of cells identified as likely for all individuals in a given sample
#' @param origin Matrix with columns indicating likely/unlikely cells for each individual.
#' @keywords assignment; stable isotopes; assignment area
#' @return Assignment area (i.e., the mean number of cells identified as likely for all individuals)
#' @examples
#' area_assign(origin = amre_assign$iso.origin, base = amre.base)
#'

area_assign <- function(origin){
  #### Average # cells identified as likely for isotope only assignment
  cells <- apply(origin,2,sum)
  area <- mean(cells)/nrow(origin)

  return(area)
}
