#' saber_forward
#'
#' Fuction repackaged from `SABER_forward_fast.R` in codebase. What's the
#'  difference with `saber_forward_parametric_conc_wise.R` ?
#'
#' @export

saber_forward <- function(
    wavelength,
    iop,
    water_type,
    view,
    sun,
    optically_shallow,
    h_w = NULL,
    rrs_bottom = NULL,
    verbose = F
  ) {

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


# RT model ----------------------------------------------------------------

  ## Remote sensing reflectance below the surface
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w
  view_w <- geometry$view_w
  rho_L <- geometry$rho_L

  ## Remote Sensing Reflectance below the water surface
  if(water_type == "1"){

    f_rs <- 0.095 # [1/sr]
  } else if (water_type == "2") {
    f_rs <- 0.0512* (1 + (4.6659 * omega_b) +
                       (-7.8387 * (omega_b^2)) +
                       (5.4571 * (omega_b^3))
    ) *
      (1 + (0.1098 / cos(sun_w))) *
      (1 + (0.4021 / cos(view_w))) # [1/sr]
  }

  rrs_0m_deep <- f_rs * omega_b # [1/sr]

  if (optically_shallow) {

    # Attenuation Coefficients
    if (water_type == 1) {
      k0 <- 1.0395 #case 1
    } else if (water_type == 2) {
      k0 <- 1.0546 #case 2
    }

    Kd <- k0*(ext/cos(sun_w))
    kuW <- (ext/cos(view_w))*((1+omega_b)^3.5421)*(1-(0.2786/cos(sun_w)))
    kuB <- (ext/cos(view_w))*((1+omega_b)^2.2658)*(1-(0.0577/cos(sun_w)))

    # Final calculation for shallow Rrs
    Ars1 <- 1.1576
    Ars2 <- 1.0389 # Parametric coeffs for shallow water

    rrs_0m_shallow <-  rrs_0m_deep *
      (1 - (Ars1 * exp(-h_w * (Kd + kuW)))) + Ars2 *
      rrs_bottom$rrs_bottom * exp(-h_w * (Kd + kuB))

    rrs_0m <- rrs_0m_shallow

  } else {
    # Optically deep water
    rrs_0m <- rrs_0m_deep

  }

  if (verbose) {
    plot(wavelength, rrs_0m, xlab="wavelength",
         ylab="Rrs 0m [m^-1]")
  }

  return(
    tibble(
      wavelength,
      rrs_0m
      )
    )

}
