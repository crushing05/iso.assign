#' Breeding assignment using weighted stable isotope and abundance data
#'
#' Function to estimate the posterior probability of origin and likely/unlikely origins for stable hydrogen isotope samples and relative breeding abundance
#' @param iso.lik Dataframe or matrix containing isotope-based likelihood of origin for each breeding cell and individual.
#' @param rel.abun  Vector of relative abundance values for each breeding cell (should sum to 1).
#' @param iso.weight  Weight value to apply to stable isotope data.
#' @param abun.weight  Weight value to apply to relative abundance data.
#' @param odds  Denominator of the odds ratio for determining likely/unlikely origins. Default is 3.
#' @keywords assignment; stable isotopes; abundance
#' @return A list containing two dataframes: prob contain posterior probabilities for each breeding origin and individual & origin contains the likely/unlikely origins for each individual
#' @examples
#' weight_assign(iso.lik = amre_iso$iso.prob, rel.abun = woth_base$rel.abun, iso.weight = -0.1, abun.weight = -0.9)
#'

abun_assign <- function(iso.like, rel.abun, iso.weight, abun.weight, odds = 3){
    iso.weight <- 10^iso.weight
    abun.weight <- 10^abun.weight
  # Estimate posteriors under each weighting combination
    prob <- prop.table(rel.abun^abun.weight*iso.like^iso.weight,2) #Prob for each indv/cell
    rel.prob <- apply(prob,2,function(x)x/max(x)) # Relative prob for each indv/cell
    origin <- ifelse(rel.prob <1/odds, 0 , 1) # Likely/unlikely for each indv/cell

  wght.summ <- list(origin = origin, prob = prob)
  return(wght.summ)
}
