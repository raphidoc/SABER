#' Lee 1998 forward model
#' Compute Rrs given initial parameters
#'
#' @author Soham Mukherjee
#'
#' @param iop a tibble with column {wavelength, a, bb}, non-water component of absorption and backscattering.
#' @param optically_shallow bolean, is the water optically shallow ? If TRUE, provide `h_w`, `rb_fraction`, `rb`.
#' @param h_w water column height [m].
#' @param rrs_botom remote sensing bottom reflectance. see function `rrs_bottom_lmm`
#' @param wavelength Requested wavelength for the simulation
#'
#' @references Lee, z. et al. (1998) ‘Hyperspectral remote sensing for shallow waters. I. A semianalytical model’, Applied Optics, 37(27), pp. 6329–6338. Available at: https://doi.org/10.1364/AO.37.006329.
#'

lee98_forward <-  function(
    wavelength,
    iop,
    view,
    sun,
    optically_shallow,
    h_w = NULL,
    rrs_bottom = NULL,
    verbose = F
){

  # Define IOPs  ------------------------------------------------------------
  iop_w <- pure_water_iop(wavelength)

  a_w <- iop_w$a_w
  bb_w <- iop_w$bb_w


  if (!identical(iop$wavelength, wavelength)) {
    rlang::abort(
      "iop walength different than requested simulation wavelength,
        see `parse_iop` or `iop_from_oac` to prepare iop data frame")

  } else {
    a_non_water <- iop$a
    bb_non_water <- iop$bb

    # Total Absorption Coefficient [1/m]
    a_t <-  a_w + a_non_water

    # Total Backscattering Coefficient [1/m]
    bb_t <-  bb_w + bb_non_water

    # extinction coeff. [1/m]
    ext <-  a_t + bb_t

    # single back scattering albedo
    omega_b <-  bb_t / ext
  }

  # Radiative Transfer Model ------------------------------------------------

  ## Remote sensing reflectance 0m the surface
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w
  view_w <- geometry$view_w
  rho_L <- geometry$rho_L

  ## Remote Sensing Reflectance 0m the water surface
  p1  <- 0.084
  p2  <- 0.17
  k1w <- 1.03
  k2w <- 2.04
  k1b <- 1.04
  k2b <- 5.04
  q1 <- 1 #(for viewing angle=0)
  k <- a_t + bb_t
  u <- bb_t / k

  rrs_0m_deep <- q1 * (p1 + p2 * u) * u

  if (optically_shallow) {

    # Attenuation Coefficients
    mu_s <- cos(sun_w)
    mu_v <- cos(view_w)
    du_w <- k1w * sqrt(1 + k2w * u)
    du_b <- k1b * sqrt(1 + k2b * u)

    rrs_0m_shallow <- rrs_0m_deep *
      (1 - exp(-(1 / mu_s + du_w / mu_v) * k * h_w)) +
      rrs_bottom$rrs_bottom * exp(-(1 / mu_s + du_b / mu_v) * k * h_w)

    rrs_0m <- rrs_0m_shallow

  } else {
    # Optically deep water
    rrs_0m <- rrs_0m_deep

  }

  if (verbose) {
    plot(wavelength, rrs_0m, xlab="wavelength",
         ylab="Rrs 0m [m^-1]")
  }

  return(tibble(
    wavelength,
    rrs_0m
  ))

}
