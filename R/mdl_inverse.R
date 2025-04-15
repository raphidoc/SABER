#' Markov Chain Monte Carlo radiative transfer inversion
#'
#' @author Soham Mukherjee
#'
#' @export

inverse_mcmc <- function(
  density,
  sampler,
  likelihood
  ) {

  # Create prior density and sampling class ---------------------------------

  if (pop.sd == "known" & type_Rrs_below == "deep") {
    prior <- BayesianTools::createPrior(
      density = density,
      sampler = sampler,
      lower = c(0,0,0),
      upper = c(30,5,0.5),
      best = as.numeric(Fit.optimized.ssobj)
      )
  }

  if (pop.sd == "unknown" & type_Rrs_below == "deep") {

    prior <- BayesianTools::createPrior(
      density = density,
      sampler = sampler,
      lower = c(0,0,0,0.0001),   # <<USER DEFINED>>
      upper = c(30,5,0.5, 0.01), # <<USER DEFINED>>
      best = c(as.numeric(Fit.optimized.ssobj),
               inverse_output[[1]]$estimates[4]))
  }

  # Create Bayesian setup for MCMC ------------------------------------------

  if (use_likelihood) {

    # With prior
    bayessetup <- BayesianTools::createBayesianSetup(
      prior = prior,
      likelihood = ll,
      #lower = c(0,0,0), upper = c(30,5,0.5),
      names = c("chl","acdom440","anap440", "pop.sd"),
      parallel = F)

  } else {

    # Only likelihood
    bayessetup <- BayesianTools::createBayesianSetup(
      prior = NULL,
      likelihood = ll,
      lower = lower.bound,
      upper = upper.bound,
      names = c("chl","acdom440","anap440", "pop.sd"),
      parallel = F)
  }

  checkBayesianSetup(bayessetup)

# Create settings for MCMC ------------------------------------------------

  settings = list(iterations = 10000, message = TRUE, nrChains = 1, burnin=2000)
  samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

  out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = samplerlist[6] )
  summary(out)

  return()
}
