input_am03 <- function(par, rrs) {
  par_requiered <- c(
    "chl",
    "a_g_440",
    "bb_p_550",
    "water_type",
    "theta_view",
    "theta_sun",
    "optically_shallow",
    "h_w")

  if (!all(par_requiered %in% names(par))
  ) {
    rlang::abort(
      glue::glue("am03 missing par: {names(par)[!par_requiered %in% names(par)]}")
    )
  }

  iop <- iop_from_oac(rrs$wavelength, par)

  r_b_fraction_vec <- par[grep("^rb_", names(par))]

  r_b <- r_b_lmm(
    r_b_spectrum = r_b_class_egsl,
    r_b_fraction = r_b_fraction_vec,
    wavelength = wavelength,
    verbose = F
  )

  list(
    wavelength = rrs$wavelength,
    iop = iop,
    water_type = par["water_type"],
    theta_view = par["theta_view"],
    theta_sun = par["theta_sun"],
    optically_shallow = par["optically_shallow"],
    h_w = par["h_w"],
    r_b = r_b
  )
}
