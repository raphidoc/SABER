#' parse_iop
#'
#' take tibbles for non water absorption and backscattering
#'
#' @author Raphael Mabit
#'
#' @param a tibble with columns {a, wavelength}
#' @param bb tibble with columns {bb, wavelength}
#' @param wavelength_out a vector with the requested interpolation wavelength
#' @param plot bollean, produce input and interpolated IOPs plot
#'
#' @returns A tibble with columns {wavelength, a, bb}
#'
#' @export

parse_iop <- function(
    a,
    bb,
    wavelength,
    verbose = F
    ) {

  a_non_water <- approx(x= a$wavelength, y = a$a,
                        xout = wavelength, method = "linear")$y

  # For bb we create a power law model and use it for interpolation

  model = nls(
    bb ~ a * wavelength ^ b,
    start = list(a = bb$bb[3], b = 1),
    data=bb,
    control = list(maxiter = 100, warnOnly = T)
    )

  # bbp_hs[j,i] = refbbp[j,1]*((waveletngth_hs[i,1]/555)^-(refexponent[j,1]))
  bb_non_water = coef(model)[1] * wavelength ^ coef(model)[2]

  # bb_non_water <- approx(x= iop$wavelength, y = iop$bb,
  #                        xout = wavelength, method = "linear")$y

  iop <- tibble(
    wavelength,
    "a" = a_non_water,
    "bb" = bb_non_water
  )

  if (verbose) {
    # Plot interpolated non water absorption
    plot(a$wavelength, a$a, xlab="wavelength",
         ylab="non-water absorption [m^-1]")
    lines(wavelength, a_non_water, col="red", lwd=3)

    # Plot power law fitted backscatering
    plot(bb$wavelength, bb$bb, xlab="wavelength",
         ylab="non-water backscatter [m^-1]")
    lines(wavelength, bb_non_water, col="red", lwd=3)
  }

  return(iop)
}
