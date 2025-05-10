generate_stan_model <- function(par_inversed, par_fixed, prior, model_code_body) {
  # Parameters block
  param_lines <- purrr::map_chr(par_inversed, function(p) {
    prior_def <- prior[[p]]
    bound_str <- glue::glue(
      "<lower={prior_def$lower}, upper={prior_def$upper}>"
    )
    glue::glue("  real{bound_str} {p};")
  })

  # Priors block
  prior_lines <- purrr::map_chr(par_inversed, function(p) {
    prior_def <- prior[[p]]
    glue::glue("  {p} ~ {prior_def$distribution};")
  })

  # Assemble complete model
  stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] rrs_obs;
  real<lower=0> sd_obs;
  // Add fixed param declarations if needed
}}

parameters {{
{paste(param_lines, collapse = '\\n')}
}}

model {{
  vector[N] rrs_model;

  // Forward model
{model_code_body}

  // Priors
{paste(prior_lines, collapse = '\\n')}

  // Likelihood
  rrs_obs ~ normal(rrs_model, sd_obs);
}}
  ")

  return(stan_code)
}

inverse_stan <- function(
    rrs,
    forward_model,
    par_inversed,
    prior = NULL,
    lower = NULL,
    best = NULL,
    upper = NULL,
    par_fixed = NULL,
    iterations = 1000,
    burnin = 500,
    chains = 4,
    stan_file = NULL
) {
  # 1. Resolve forward model and STAN file
  model_entry <- get_forward_model(forward_model)
  if (is.null(stan_file)) {
    stan_file <- model_entry$stan_file
  }

  # 2. Combine fixed + variable parameters
  par_all <- c(par_inversed, names(par_fixed %||% list()))

  # 3. Prepare STAN data
  stan_data <- model_entry$prepare_stan_data(list(
    rrs_observed = rrs,
    par_fixed = par_fixed
  ))

  # 4. Compile and run STAN
  mod <- cmdstanr::cmdstan_model(
    stan_file = "stan/am03_model.stan",
    user_header = "src/am03_model.hpp",
    cpp_options = list(stan_threads = TRUE)
  )

  fit <- mod$sample(
    data = stan_data,
    iter_warmup = burnin,
    iter_sampling = iterations,
    chains = chains,
    parallel_chains = chains
  )

  # 5. Extract posterior summary
  draws <- fit$draws_df()
  par_mean <- posterior::summarise_draws(draws)[, c("variable", "mean", "sd")]
  par_mean <- dplyr::filter(par_mean, variable %in% par_inversed)

  estimates <- stats::setNames(
    c(par_mean$mean, par_mean$sd),
    c(par_mean$variable, paste0(par_mean$variable, "_sd"))
  )

  return(estimates)
}
