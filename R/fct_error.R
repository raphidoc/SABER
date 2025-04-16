#' log likelihood
#'
#' @return the error (or likelihood) of modelled vs observed values

log_ll <- function(modelled, observed, sd) {
  error <- -sum(
    dnorm(x = 10000 * observed, mean = 10000 * modelled, sd = sd, log = TRUE)
  )
  return(error)
}

#' residual sum of square
#'
#' @return the error (or likelihood) of modelled vs observed values

rss <- function(modelled, observed) {
  sum((observed - modelled)^2, na.rm = TRUE)
}

#' lee 1998 spectral error index
#'
#' @return the error (or likelihood) of modelled vs observed values
#' @references Lee, Z. et al. (1999) ‘Hyperspectral remote sensing for shallow waters: 2 Deriving bottom depths and water properties by optimization’, Applied Optics, 38(18), p. 3831. Available at: https://doi.org/10.1364/AO.38.003831.

lee99 <- function(modelled, observed, wavelength) {
  region1 <- which(wavelength >= 400 & wavelength <= 675)
  region2 <- which(wavelength >= 750 & wavelength <= 830)
  numerator <- sqrt(sum((observed[region1] - modelled[region1])^2) +
    sum((observed[region2] - modelled[region2])^2))
  denominator <- sum(modelled[region1]) + sum(modelled[region2])
  error <- numerator / denominator
  return(error)
}
