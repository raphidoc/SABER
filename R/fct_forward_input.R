input_am03 <- function(par, rrs) {
  # par_requiered <- c(
  #   "chl",
  #   "a_g_440",
  #   "bb_p_550",
  #   "water_type",
  #   "theta_view",
  #   "theta_sun",
  #   "optically_shallow",
  #   "h_w")
  #
  # if (!all(par_requiered %in% names(par))
  # ) {
  #   rlang::abort(
  #     glue::glue("am03 missing par: {names(par)[!par_requiered %in% names(par)]}")
  #   )
  # }

  iop <- iop_from_oac(rrs$wavelength, par)

  r_b_fraction_vec <- par[grep("^r_rs_b", names(par))]

  r_b <- compute_r_rs_b_lmm(
    fractions = r_b_fraction_vec
  )

  list(
    wavelength = rrs$wavelength,
    iop = iop,
    water_type = par["water_type"],
    theta_view = par["theta_view"],
    theta_sun = par["theta_sun"],
    h_w = par["h_w"],
    r_b = r_b
  )
}
