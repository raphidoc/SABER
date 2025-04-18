make_prior <- function(sample_data, dist = "weibull", name = "param") {
  fit <- fitdistrplus::fitdist(sample_data, dist)
  list(
    fit = fit,
    density = function(x) {
      do.call(paste0("d", dist), list(x = x, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE))
    },
    sampler = function(n = 1) {
      do.call(paste0("r", dist), list(n = n, shape = fit$estimate[1], scale = fit$estimate[2]))
    },
    name = name
  )
}


make_prior_bundle <- function(priors, lower, upper, best_guess = NULL) {
  density_fn <- function(par) {
    sum(mapply(function(p, val) p$density(val), priors, par))
  }

  sampler_fn <- function(n = 1) {
    matrix(unlist(lapply(priors, function(p) p$sampler(n))), ncol = length(priors))
  }

  list(
    prior = BayesianTools::createPrior(
      density = density_fn,
      sampler = sampler_fn,
      lower = lower,
      upper = upper,
      best = best_guess
    ),
    names = sapply(priors, function(p) p$name)
  )
}

inverse_mcmc <- function(
    prior_bundle,
    likelihood,
    iterations = 10000,
    burnin = 2000,
    sampler = "DEzs",
    use_prior = TRUE
    ) {

  # Wrap likelihood to restore names
  prior_names <- prior_bundle$names

  named_likelihood <- function(par) {
    names(par) <- prior_names
    likelihood(par)
  }

  setup <- BayesianTools::createBayesianSetup(
    prior = if (use_prior) prior_bundle$prior else NULL,
    likelihood = named_likelihood,
    lower = c(0, 0, 0),
    upper = c(100, 5, 0.1),
    names = prior_names,
    parallel = FALSE
  )

  BayesianTools::checkBayesianSetup(setup)

  out <- BayesianTools::runMCMC(
    bayesianSetup = setup,
    settings = list(
      iterations = iterations,
      burnin = burnin,
      message = TRUE,
      nrChains = 1
    ),
    sampler = sampler
  )

  return(out)
}


#' #' density_fct
#' #'
#' #' density function input pareter to `BayesianTools::createPrior`
#' #'
#'
#' density <- function(par) {
#'   chl <- par[1]
#'   acdom.440 <- par[2]
#'   anap.440 <- par[3]
#'   x.sd <- par[5]
#'
#'   density_chl <- dweibull(
#'     x = chl,
#'     shape = fit.chl.norm$estimate[1],
#'     scale = fit.chl.norm$estimate[2],
#'     log = T
#'   )
#'
#'   density_ag440 <- dweibull(
#'     x = acdom.440,
#'     shape = fit.acdom440.norm$estimate[1],
#'     scale = fit.acdom440.norm$estimate[2],
#'     log = T
#'   )
#'
#'   density_anap440 <- dweibull(
#'     x = anap.440,
#'     shape = fit.anap440.norm$estimate[1],
#'     scale = fit.anap440.norm$estimate[2],
#'     log = T
#'   )
#'
#'   density_lklhood <- dunif(
#'     x = x.sd,
#'     min = 0.0001,
#'     max = 0.01,
#'     log = T
#'   )
#'
#'   return(density_chl + density_ag440 + density_anap440 + density_lklhood)
#' }
#'
#'
#' #' sampler_fct
#' #'
#' #' sampler function input pareter to `BayesianTools::createPrior`
#' #'
#'
#' sampler_fct <- function(n = 1) {
#'   density_chl <- rweibull(
#'     n,
#'     shape = fit.chl.norm$estimate[1],
#'     scale = fit.chl.norm$estimate[2]
#'   )
#'
#'   density_ag440 <- rweibull(
#'     n,
#'     shape = fit.acdom440.norm$estimate[1],
#'     scale = fit.acdom440.norm$estimate[2]
#'   )
#'
#'   density_anap440 <- rweibull(
#'     n,
#'     shape = fit.anap440.norm$estimate[1],
#'     scale = fit.anap440.norm$estimate[2]
#'   )
#'
#'   density_lklhood <- runif(
#'     n,
#'     min = 0.0001,
#'     max = 0.01
#'   )
#'
#'   return(
#'     cbind(
#'       density_chl,
#'       density_ag440,
#'       density_anap440,
#'       density_lklhood
#'     )
#'   )
#' }
