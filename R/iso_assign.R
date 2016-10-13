#' Breeding assignment using stable isotope data
#'
#' Function to estimate the likelihood of origin and likely/unlikely origins for stable hydrogen isotope samples
#' @param dd  Vector of length = # individuals containing the feather deuterium values.
#' @param df_base  Vector containing the predicted deuterium values for each potential breeding location.
#' @param lat Vector containing the latitudes for each potential breeding location.
#' @param lon Vector containing the longitudes for each potential breeding location.
#' @param odds  Odds for determining likely/unlikely origins. Default is 0.67.
#' @param names Optional vector containing individual IDs
#' @keywords assignment; stable isotopes
#' @return A data frame containing:
#'  indv: Individual ID
#'  lat: Latitude of each breeding cell
#'  lon: Longitude of each breeding cell
#'  iso_cut: The individual probabiilty threshold for defining likely/unlikely cells
#'  iso.like: The likelihood values for each breeding cell
#'  iso.prob: The relative probabilities for each breeding cell
#'  iso.origin: The likely/unlikely origins for each individual
#' @examples
#' iso_assign(dd = woth_dd, df_base = woth_base$ddf, lat = woth_base$y, lon = woth_base$x)
#'

iso_assign <- function(dd, df_base, lat, lon, odds = 0.67, names) {
  #### Create empty matrix to store likelihoods
  iso_like <- matrix(ncol = length(dd), nrow = length(df_base)) # One row/cell; One column/indv
  if(missing(names)){
    colnames(iso_like) <- paste("Indv_", seq(1:length(dd)), sep = "")
  }else{
    colnames(iso_like) <- names
  }

  #### Estimate likelihoods
  for(i in 1:length(dd)){ # For each individual
    iso_like[, i] <- dnorm(dd[i], mean = df_base, sd = 12)
  }


  #### Create output object
  iso_data <- as.data.frame(iso_like) %>%       # Convert likelihood matrix to data frame
    ## Covert to long format
    tidyr::gather(indv, iso_like) %>%
    ## For each individual, convert likelihood's to probabilities (sum to 1), add basemap lat/long
    dplyr::group_by(indv) %>%
    dplyr::mutate(iso_prob = iso_like / sum(iso_like),
                  lat = lat, lon = lon) %>%
    dplyr::ungroup() %>%
    ## Nest by individual
    tidyr::nest(-indv) %>%
    ## Estimate cumulative probability and estimate cutoff value for each individual
    dplyr::mutate(fit = purrr::map(data, ~ predict(smooth.spline(cumsum(sort(.$iso_prob)), sort(.$iso_prob), spar = 0.1), 1 - odds)),
                  iso_cut = map_dbl(fit, "y")) %>%
    ## Remove spline predictions & unnest dataframe
    dplyr::select(-fit) %>%
    tidyr::unnest() %>%
    ## Reclassify cells as likely/unlikely based on cumulative prob
    dplyr::mutate(iso_origin = ifelse(iso_prob > iso_cut, 1, 0)) %>%
    ## Reorder columns
    dplyr::select(indv, lat, lon, iso_cut, iso_like, iso_prob, iso_origin)

  return(iso_data)
}
