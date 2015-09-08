#' Breeding assignment using stable isotope data
#'
#' Function to estimate the likelihood of origin and likely/unlikely origins for stable hydrogen isotope samples
#' @param dd  Vector of length = # individuals containing the feather deuterium values.
#' @param df.base  Vector of predicted deuterium values for each potential breeding location.
#' @param odds  Denominator of the odds ratio for determining likely/unlikely origins. Default is 3.
#' @keywords assignment; stable isotopes
#' @return A list containing three dataframes:
#'  iso.like contains a likelihood values for each breeding origin and each individual
#'  iso.prob contain relative probabilities for each breeding origin and individual
#'  iso.origin contains the likely/unlikely origins for each individual
#' @examples
#' iso_assign(dd = woth_dd, df.base = woth_base$ddf)
#'

iso_assign <- function(dd, df.base, odds = 3) {
  #### Create empty data.frame to store likelihoods
  assign <- data.frame(matrix(ncol=length(dd),nrow=length(df.base))) # One row/cell; One column/indv
  colnames(assign) <- seq(1:length(dd))
  tassign <- t(assign)
  colnames(tassign) <- round(df.base,digits = 2)
  assign <- t(tassign)

  #### Estimate likelihoods
    for(i in 1:length(dd)){ # For each individual
      for(k in 1:length(df.base)){  # For each basemap cell
        assign[k,i]=dnorm(dd[i], mean = df.base[k], sd = 12)
      }
    }

  iso.prob <- apply(assign,2,function(x) x/max(x)) # Covert to relative probabilities
  iso.origin <- ifelse(iso.prob < 1/odds, 0, 1)    # Convert to likely/unlikely origin


  #### Store data as list and return
  iso.data <- list(iso.like = data.frame(assign), iso.prob = data.frame(iso.prob), iso.origin = data.frame(iso.origin))
  return(iso.data)
}
