log_ll <- function(modelled, observed, sd) {
  error <- -sum(
    dnorm(x = 10000 * observed, mean = 10000 * modelled, sd = sd, log = TRUE)
    )

  return(error)
}

lee99 <- function(modelled, observed, wavelength) {
  region1 <- which(wavelength >= 400 & wavelength <= 675)
  region2 <- which(wavelength >= 750 & wavelength <= 830)
  numerator <- sqrt(sum((observed[region1] - modelled[region1])^2) +
                      sum((observed[region2] - modelled[region2])^2))
  denominator <- sum(modelled[region1]) + sum(modelled[region2])
  error <- numerator / denominator

  return(error)
}
