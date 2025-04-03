#' saber_inverse
#' Implement the SABER inversion. BGC stands for biogeochemstry
#'
#' @param sa_model c("am03", "lee98")
#' @param objective_fct c("log-LL", "SSR", "obj_L98")
#' @param method_opt c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent","levenberg-marqardt", "auglag")
#'

solve_objective_inverse_shallow <- function(
    constrain_shallow_iop,
    abs_path,
    bb_path,
    constrained_bgc,
    constrain_bgc_param = c("chl", "adg", "bbp"),
    constrain_bgc_value,
    unconstrained,
    auto_spectral_slope,
    manual_spectral_slope,
    manual_spectral_slope_vals = c("s_g"=0.015, "s_d"=0.01160, "gamma"=0.5),
    initial,
    objective_fct = c("log-LL", "SSR", "obj_L98"),
    sa_model="am03",
    method_opt = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                   "Brent","levenberg-marqardt", "auglag"),
    obsdata,
    wavelength,
    lower_b,
    upper_b
){

  initial_rb_length = 3


  if (auto_spectral_slope == TRUE & manual_spectral_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  myFun <- function(x) {
    NA
  }








}
