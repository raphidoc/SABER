#' fct_compose_minimization
#' a function to dynamically create the minimization function
#'  passed to inversion
#'
#' @author Raphael Mabit
#' @param forward_model c("am03", "lee98")
#' @param error_fct c("log-ll", "lee99")

fct_compose_minimization <- function(
    forward_model,
    error_fct,
    rrs_observed,
    par_fixed) {
  # TODO: should be a packages variable lazy loaded
  forward_models <- c("am03", "lee98")
  error_fcts <- c("log-ll", "lee99")


  # Select forward model ----------------------------------------------------
  forward_model <- switch(forward_model,
    "am03" = saber_forward,
    "lee98" = lee98_forward,
    rlang::abort(
      glue::glue("forward_model must be one of:", forward_models)
    )
  )


  # Prepare generic model inputs --------------------------------------------
  # Internal function to combine dynamic and fixed parameters
  complete_par <- function(par) {
    if (!is.null(par_fixed)) {
      # Ensure no overlap
      if (any(par_fixed$parameter %in% names(par))) {
        stop("Parameters cannot be defined in both `par_init` and `par_fixed`.")
      }

      fixed_vals <- stats::setNames(par_fixed$value, par_fixed$parameter)
      par <- c(par, fixed_vals)
    }
    par[order(names(par))]
  }

  prepare_model_inputs <- function(par) {
    oac <- tibble(
      chl = par["chl"],
      ag_440 = par["ag_440"],
      bbp_550 = par["bbp_550"]
    )

    iop <- iop_from_oac(
      rrs_observed$wavelength,
      oac
    )

    # extract parameter starting with "rb_" (fractions of rrs_end_member)
    end_member_fraction <- par[grep("^rb_", names(par))]

    rrs_bottom <- rrs_bottom_lmm(
      rrs_end_member = rb_algae_wise,
      end_member_fraction = end_member_fraction,
      wavelength = rrs_observed$wavelength,
      verbose = FALSE
    )

    return(list(iop = iop, rrs_bottom = rrs_bottom))
  }


  # Construct minimization function -----------------------------------------
  minimization_fct <- switch(error_fct,
    "log-ll" = function(par) {
      par <- complete_par(par)
      inputs <- prepare_model_inputs(par)

      rrs_modeled <- forward_model(
        wavelength = rrs_observed$wavelength,
        iop = inputs$iop,
        view = 0,
        sun = 20,
        optically_shallow = T,
        h_w = par["h_w"],
        rrs_bottom = inputs$rrs_bottom,
        verbose = F
      )

      return(log_ll(
        modelled = rrs_modeled$rrs_0m,
        observed = rrs_observed$rrs_0m,
        sd = par["sd"]
      ))
    },
    "lee99" = function(par) {
      par <- complete_par(par)
      inputs <- prepare_model_inputs(par)

      rrs_modeled <- forward_model(
        wavelength = rrs_observed$wavelength,
        iop = inputs$iop,
        view = 0,
        sun = 20,
        optically_shallow = T,
        h_w = par["h_w"],
        rrs_bottom = inputs$rrs_bottom,
        verbose = F
      )

      return(lee99(
        modelled = rrs_modeled$rrs_0m,
        observed = rrs_observed$rrs_0m,
        wavelength = rrs_observed$wavelength
      ))
    }
  )

  return(minimization_fct)
}
