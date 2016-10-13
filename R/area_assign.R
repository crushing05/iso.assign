#' Estimate the assignment area for known-origin tissue samples
#'
#' Function to estimate the mean number of cells identified as likely for all individuals in a given sample
#' @param summ Assignment results from either the iso_assign() or abun_assign() functions
#' @param iso Should area be estimated for the isotope-only assignment (default) or with abundance used a prior
#' @keywords assignment; stable isotopes; assignment area
#' @return Assignment area (i.e., the mean number of cells identified as likely for all individuals)
#' @examples
#' area_assign(summ = woth_iso)
#'

area_assign <- function(summ, iso = TRUE){
  #### Average # cells identified as likely for isotope only assignment
  if(iso == TRUE){
    area <- summ %>%
      group_by(indv) %>%
      summarise(cells = sum(iso_origin), area = cells / length(iso_origin)) %>%
      .$area
  }else{
    area <- summ %>%
      group_by(indv) %>%
      summarise(cells = sum(wght_origin), area = cells / length(wght_origin)) %>%
      .$area
  }

  return(mean(area))
}
