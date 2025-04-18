.onLoad <- function(libname, pkgname) {
  register_input_preparer("input_am03", function(par, rrs) {
    input_am03(par, rrs)
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
      r_b = par$r_b
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
      r_b = par$r_b
    )
  })

  register_objective_function("log-ll", function(modelled, observed, par) {
    log_ll(modelled = modelled, observed = observed, sd = par[["sd"]])
  })

  register_objective_function("rss", function(modelled, observed, par) {
    rss(modelled = modelled, observed = observed)
  })

  register_objective_function("lee99", function(modelled, observed, par) {
    lee99(modelled = modelled, observed = observed, wavelength = par[["wavelength"]])
  })
}
