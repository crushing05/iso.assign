#' Breeding assignment using weighted stable isotope and abundance data
#'
#' Function to estimate the posterior probability of origin and likely/unlikely origins for stable hydrogen isotope samples and relative breeding abundance
#' @param summ Dataframe containing the isotope- and abundance-base assignment results from the iso_assign() and/or abun_assign() function
#' @param iso  Should coordinates be based on isotope assignment or abundance assignment?
#' @keywords assignment; stable isotopes; abundance
#' @return Data frame containing the estimate latitude & logitude and uncertainty for each individual.
#' @examples
#' woth_coord <- wght_coord(summ = woth_iso)
#'

### Computed mean lat/long, weighted by posteriors

wght_coord <- function(summ, iso = TRUE) {
  if(iso == TRUE){
    summ <- rename(summ, origin = iso_origin, prob = iso_prob)
  }else{
    summ <- rename(summ, origin = wght_origin, prob = wght_prob)
  }

  coords <- summ %>%
    mutate(prob = prob * origin, lat = lat * origin, lon = lon * origin,
           lat1 = lat * pi / 180, lon1 = lon * pi / 180,
           X1 = cos(lat1) * cos(lon1), Y1 = cos(lat1) * sin(lon1), Z1 = sin(lat1)) %>%
    group_by(indv) %>%
    summarise(nCell = sum(origin),
              V1 = sum(prob),
              V2 = sum(prob ^ 2),
              varX = sum((X1[X1 != 1] - mean(X1[X1 != 1])) ^ 2),
              varY = sum((Y1[Y1 != 0] - mean(Y1[Y1 != 0])) ^ 2),
              varZ = sum((Z1[Z1 != 0] - mean(Z1[Z1 != 0])) ^ 2),
              wght_x = sum(prob * X1) / V1,
              wght_y = sum(prob * Y1) / V1,
              wght_z = sum(prob * Z1) / V1,
              var_wght_x = (1 - V2 / V1) * varX,
              var_wght_y = (1 - V2 / V1) * varY,
              var_wght_z = (1 - V2 / V1) * varZ,
              wght_x_l = wght_x - 1.96 * sqrt(var_wght_x / nCell),
              wght_y_l = wght_y - 1.96 * sqrt(var_wght_y / nCell),
              wght_z_l = wght_z - 1.96 * sqrt(var_wght_z / nCell),

              wght_x_u = wght_x + 1.96 * sqrt(var_wght_x / nCell),
              wght_y_u = wght_y + 1.96 * sqrt(var_wght_y / nCell),
              wght_z_u = wght_z + 1.96 * sqrt(var_wght_z / nCell),


              # Convert average coordinates to lat/long
              Lon = atan2(x = wght_x, y = wght_y),
              Lat = asin(wght_z),

              Lon_l = atan2(x = wght_x_l, y = wght_y_l),
              Lat_l = asin(wght_z_l),

              Lon_u = atan2(x = wght_x_u, y = wght_y_u),
              Lat_u = asin(wght_z_u),


              # Convert lat/long to degrees
              x = Lon * 180 / pi,
              y = Lat * 180 / pi,

              x_l = Lon_l * 180 / pi,
              y_l = Lat_l * 180 / pi,

              x_u = Lon_u * 180 / pi,
              y_u = Lat_u * 180 / pi) %>%
    select(indv, x, y, x_l, y_l, x_u, y_u) %>%
    rename(lon = x, lat = y, lon_LCI = x_l, lat_LCI = y_l, lon_UCI = x_u, lat_UCI = y_u)
  return(coords)
}

