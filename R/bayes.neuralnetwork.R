library(brms)

# Generate some simulated data
set.seed(123)

#Load NOMAD in situ Prior data
bgc.data <- read.csv("./nomad_dataset_simplified.csv",
                     header = T, sep = ",")


chl.sample <- bgc.data$chl
chl.sample = chl.sample[!chl.sample %in% -999]

#hist(chl.sample, probability = T)

acdom.440.sample <- bgc.data$ag443
acdom.440.sample = acdom.440.sample[!acdom.440.sample %in% -999]

#hist(acdom.440.sample, probability = T)

anap.440.sample <- bgc.data$ad443
anap.440.sample = anap.440.sample[!anap.440.sample %in% -999]

#hist(anap.440.sample, probability = T)

fit.chl.norm <- fitdistrplus::fitdist(chl.sample, "weibull") #try to fit dist

# curve(dweibull(x, shape = fit.chl.norm$estimate["shape"], 
#                                scale = fit.chl.norm$estimate["scale"] ), add = T)

#plot(fit.chl.norm) #Infer results of the fit

fit.acdom440.norm <- fitdistrplus::fitdist(acdom.440.sample, "weibull") #try to fit dist

#plot(fit.acdom440.norm) #Infer results of the fit

fit.anap440.norm <- fitdistrplus::fitdist(anap.440.sample, "weibull") #try to fit dist
#plot(fit.anap440.norm) #Infer results of the fit

# Define the likelihood function
likelihood <- function(chl, acdom440, anap440, pop_sd) {
  Gpred = Saber_forward(chl = chl, acdom440 = acdom440,
                        anap440 = anap440,
                        #bbp.550 = HL.deep.iop$bbp550[j],
                        bbp.550 = Fit.input$bbp.550,
                        verbose=F ,realdata = obsdata)
  

  
  smull = dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pop_sd, 
                    log = TRUE)  #10000 is the scaling factor to avoid calculation  
                                 #of very small numbers
  
  
  return(smull)
}

likelihood(chl = 2, acdom440 = 2, anap440 = 1, pop_sd=0.001)

# Define the prior distributions
prior_chl <- brms::prior(prior = 'weibull', shape = as.numeric(fit.chl.norm$estimate[1]), 
                         scale =as.numeric(fit.chl.norm$estimate[2]))

prior_acdom440 <- brms::prior('weibull', shape = fit.acdom440.norm$estimate[1], 
                              scale =fit.acdom440.norm$estimate[2], nl=n_obs)

prior_anap440 <- brms::prior('weibull', shape = fit.anap440.norm$estimate[1], 
                             scale =fit.anap440.norm$estimate[2], nl=n_obs)

#examine the prior options and the brms default
prior = get_prior(habitat_lost ~ mass + diet + (1|IUCN), data = data, 
                  family = weibull())

prior_chl <- prior('weibull', shape=1.5, scale=1, nl=n_obs)
prior_acdom440 <- prior('weibull', shape=1.5, scale=1, nl=n_obs)
prior_anap440 <- prior('weibull', shape=1.5, scale=1, nl=n_obs)

# Define the function for the SA model
sa_model <- function(chl, acdom440, anap440) {
  # Define the SA model of Albert & Mobley 2003 here
  return(rrs)
}

# Define the neural network to approximate the SA model
nn_model <- function(input_dim, hidden_dim, output_dim) {
  brm(
    formula = bf(y ~ nn(hidden(dim=input_dim, activation='tanh') | chl + acdom440 + anap440)),
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), nl=hidden_dim, coefname=paste0('nn_', 1:hidden_dim, '_hidden_')),
      prior(normal(0, 1), nl=output_dim, coefname=paste0('nn_hidden_', 1:output_dim, '_output'))
    ),
    chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95, max_treedepth = 15),
    data = list(chl = chl_true, acdom440 = acdom440_true, anap440 = anap440_true, y = y)
  )
}

# Define the Bayesian Deep Learning model
bdl_model <- function(input_dim, hidden_dim, output_dim) {
  brm_multiple(
    formula = likelihood,
    family = gaussian(),
    prior = c(prior_chl, prior_acdom440, prior_anap440),
    data = list(y = y),
    nn = list(nn_model(input_dim, hidden_dim, output_dim))
  )
}

# Run the MCMC algorithm to estimate the posterior distribution
set.seed(123)
bdl <- bdl_model(input_dim=3, hidden_dim=10, output_dim=1)
posterior <- posterior_samples(bdl, pars = c('chl', 'acdom440', 'anap440', 'nn_y'))

# Use the posterior distribution to make predictions with uncertainty estimates
chl_post <- posterior$chl
acdom440_post <- posterior$acdom440
anap440_post <- posterior$anap440
y_post <- posterior$nn_y
