#' Make likelihood/minimization function
#'
#' @param model Name of the forward model
#' @param error Name of the error function
#' @param rrs_observed Observed Rrs data
#' @param par_fixed Named list or data frame of fixed parameters
#'
#' @return A function returning the model error/likelihood
#' @export
make_likelihood <- function(model, error, rrs_observed, par_fixed = NULL) {
  prepare_input <- get_input_preparer(paste0("input_", model))
  forward_model <- get_forward_model(model)
  error_function <- get_error_function(error)

  if (is.null(prepare_input)) stop(paste("Unknown prepare input for: ", paste0("input_", model)))
  if (is.null(forward_model)) stop(paste("Unknown forward model: ", model))
  if (is.null(error_function)) stop(paste("Unknown error function: ", error))

  complete_par <- function(par) {
    if (!is.null(par_fixed)) {

      if (any(par_fixed$parameter %in% names(par))) {
        stop("Parameters cannot be defined in both `par_init` and `par_fixed`.")
      }

      fixed_vals <- stats::setNames(par_fixed$value, par_fixed$parameter)
      par <- c(par, fixed_vals)
    }

    par[order(names(par))]
  }

  function(par) {
    par <- complete_par(par)
    inputs <- prepare_input(par, rrs_observed)
    rrs_modeled <- forward_model(inputs)

    error_function(
      modelled = rrs_modeled$rrs_0m,
      observed = rrs_observed$rrs_0m,
      par = par
    )
  }
}
