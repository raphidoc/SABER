#' @title Forward model input preparation Registry
#' @description List of forward modelsinput preparation available to SABER
.input_preparer_registry <- new.env(parent = emptyenv())

#' @title Forward Model Registry
#' @description List of forward models available to SABER
.forward_model_registry <- new.env(parent = emptyenv())

#' @title Error Function Registry
#' @description List of error functions available to SABER
.error_function_registry <- new.env(parent = emptyenv())

#' Register a new input preparation
#' @param name Name of the input preparation (must end with *input_\*modelName\**)
#' @param fn A function taking model parameters and returning model inputs
register_input_preparer <- function(name, fn) {
  .input_preparer_registry[[name]] <- fn
}

#' Register a new forward model
#' @param name Name of the model
#' @param fn A function taking model parameters and returning modeled Rrs
register_forward_model <- function(name, fn) {
  .forward_model_registry[[name]] <- fn
}

#' Register a new error function
#' @param name Name of the error function
#' @param fn A function that returns a single numeric error value
register_error_function <- function(name, fn) {
  .error_function_registry[[name]] <- fn
}

#' Get registered model or error function
get_input_preparer <- function(name) {
  get0(name, envir = .input_preparer_registry)
}

get_forward_model <- function(name) {
  get0(name, envir = .forward_model_registry)
}

get_error_function <- function(name) {
  get0(name, envir = .error_function_registry)
}

#' List all registered input preparer
#' @export
list_input_preparer <- function() {
  ls(.input_preparer_registry)
}

#' List all registered forward models
#' @export
list_forward_models <- function() {
  ls(.forward_model_registry)
}

#' List all registered error functions
#' @export
list_error_functions <- function() {
  ls(.error_function_registry)
}
