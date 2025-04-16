.onLoad <- function(libname, pkgname) {
  register_input_preparer("input_am03", function(par, rrs) {

    par_requiered <- c(
      "chl",
      "ag_440",
      "bbp_550",
      "water_type",
      "theta_view",
      "theta_sun",
      "optically_shallow",
      "h_w")

    if (!all(par_requiered %in% names(par))
    ) {
      rlang::abort(
        glue::glue("Missing par: {names(par)[!par_requiered %in% names(par)]}")
        )
    }

    oac <- tibble::tibble(
      chl = par["chl"],
      ag_440 = par["ag_440"],
      bbp_550 = par["bbp_550"]
    )

    iop <- iop_from_oac(rrs$wavelength, oac, rrs)

    end_member_fraction <- par[grep("^rb_", names(par))]
    # TODO: better manage rrs_end_member
    rrs_bottom <- rrs_bottom_lmm(
      rrs_end_member = rb_algae_wise,
      end_member_fraction = end_member_fraction,
      wavelength = rrs$wavelength,
      verbose = FALSE
    )

    list(
      wavelength = rrs$wavelength,
      iop = iop,
      water_type = par["water_type"],
      theta_view = par["theta_view"],
      theta_sun = par["theta_sun"],
      optically_shallow = par["optically_shallow"],
      h_w = par["h_w"],
      rrs_bottom = rrs_bottom
    )
  })


  register_forward_model("am03", function(par) {
    forward_am03(
      wavelength = par$wavelength,
      iop = par$iop,
      water_type = par$water_type,
      theta_view = par$theta_view,
      theta_sun = par$theta_sun,
      optically_shallow = par$optically_shallow,
      h_w = par$h_w,
      rrs_bottom = par$rrs_bottom
    )
  })

  register_forward_model("lee98", function(par) {
    forward_lee98(
      wavelength = par$wavelength,
      iop = par$iop,
      theta_view = par$theta_view,
      theta_sun = par$theta_sun,
      optically_shallow = par$optically_shallow,
      h_w = par$h_w,
      rrs_bottom = par$rrs_bottom
    )
  })

  register_error_function("log-ll", function(modelled, observed, par) {
    log_ll(modelled = modelled, observed = observed, sd = par[["sd"]])
  })

  register_error_function("rss", function(modelled, observed, par) {
    rss(modelled = modelled, observed = observed)
  })

  register_error_function("lee99", function(modelled, observed, par) {
    lee99(modelled = modelled, observed = observed, wavelength = par[["wavelength"]])
  })
}
