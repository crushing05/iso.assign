#' Breeding assignment using weighted stable isotope and abundance data
#'
#' Function to estimate the posterior probability of origin and likely/unlikely origins for stable hydrogen isotope samples and relative breeding abundance
#' @param iso_data Dataframe containing the isotope-based assignment results from the iso_assign() function
#' @param rel_abun  Vector of relative abundance values for each breeding cell (should sum to 1).
#' @param iso_weight  Weight value to apply to stable isotope data.
#' @param abun_weight  Weight value to apply to relative abundance data.
#' @param odds  Odds for determining likely/unlikely origins. Default is 0.67.
#' @keywords assignment; stable isotopes; abundance
#' @return Data frames containing the isotope-based assignment results plus the weighted posterior probabilities & origins for each individual
#' @examples
#' woth_iso <- iso_assign(dd = woth_dd, df_base = woth_base$ddf, lat = woth_base$y, lon = woth_base$x)
#' abun_assign(iso_Data = woth_iso, rel_abun = woth_base$rel.abun, iso_weight = -0.7, abun_weight = 0)
#'

abun_assign <- function(iso_data, rel_abun, iso_weight, abun_weight, odds = 0.67){
    iso.weight <- 10^iso_weight
    abun.weight <- 10^abun_weight

  # Estimate posteriors under each weighting combination
    wght_summ <- iso_data %>%
      dplyr::group_by(indv) %>%
      dplyr::mutate(temp_prob = rel_abun ^ abun.weight * iso_like ^ iso.weight,
                    wght_prob = temp_prob / sum(temp_prob)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-temp_prob) %>%
      ## Nest by individual
      tidyr::nest(-indv) %>%
      ## Estimate cumulative probability and estimate cutoff value for each individual
      dplyr::mutate(fit = purrr::map(data, ~ predict(smooth.spline(cumsum(sort(.$wght_prob)), sort(.$wght_prob), spar = 0.1), 1 - odds)),
                    wght_cutoff = map_dbl(fit, "y")) %>%
      ## Remove spline predictions & unnest dataframe
      dplyr::select(-fit) %>%
      tidyr::unnest() %>%
      ## Reclassify cells as likely/unlikely based on cumulative prob
      dplyr::mutate(wght_origin = ifelse(wght_prob > wght_cutoff, 1, 0))

  return(wght_summ)
}
