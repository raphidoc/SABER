#=============================================================================================
# Functions necessary to implement MCMC sampling of parameter of interest to retrieve the water 
#components along with bathymetry and bottom reflectance (for shallow water) given the 
#parameter range and Rrs as only inputs.

#The functions developed below are following :

#1. Function to fit different distribution to the parameters of interest, i.e. Phi(par).
#2. Function to model the prior [P(theta)] using the optimal distribution(s) of Phi(par). 
#   [both desity and random sampler]
#3. Function to model the conditional log-likelihood (P(Data|Theta) of the model noise.
#4. Estimate Posterior; P(theta|Data) ~= log(P(theta)) + log(P(Data|Theta)

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================================

#--------------------------------------------------------------------------
# Function to fit different distributions to observed data for deep waters
#--------------------------------------------------------------------------
create.prior.data <- function(use.ioccg.prior=F, use.wise.prior=F, use.nomad.prior=T,
                              plot.diag=F, 
                              distrib.fit = c("weibull", "lognormal", "normal"),
                              truncate_chl = c(0.5,30),
                              truncate_adg = c(0.1,5),
                              truncate_bbp = c(0.002,0.01),
                              sample_count = 100
                              ){
  
  if (use.ioccg.prior == TRUE) {
    
    HL.deep.iop = read_IOCCG_data()
    chl.sample = HL.deep.iop$chl
    adg.sample = HL.deep.iop$acdom440+HL.deep.iop$anap440
    bbp.sample = HL.deep.iop$bbp550
    
    #Truncate the vector of prior parameters
    if (!all(is.na(truncate_chl))) {
      chl.sample = chl.sample[chl.sample >= truncate_chl[1] & chl.sample <= truncate_chl[2]]
      
    }
    
    if (!all(is.na(truncate_adg))) {
      adg.sample = adg.sample[adg.sample >= truncate_adg[1] & adg.sample <= truncate_adg[2]]
      
    }
    
    if (!all(is.na(truncate_bbp))) {
      bbp.sample = bbp.sample[bbp.sample >= truncate_bbp[1] & bbp.sample <= truncate_bbp[2]]
      
    }
    
    if (!is.na(sample_count)) {
      set.seed(123)
      chl.sample = sample(chl.sample, size = sample_count)
      adg.sample = sample(adg.sample, size = sample_count)
      bbp.sample = sample(bbp.sample, size = sample_count)
    }
    
    #try to fit the user-defined distribution
    fit.chl.norm <- fitdistrplus::fitdist(chl.sample, distrib.fit) 
    
    fit.acdm440.norm <- fitdistrplus::fitdist(adg.sample,distrib.fit) 
    
    fit.bbp550.norm <- fitdistrplus::fitdist(bbp.sample,distrib.fit) 
    
  } 
  
  if (use.wise.prior == TRUE) {
    
    #Load WISE-Man 2019 in situ Prior data
    bgc.data <- read.csv("./data/biogeochemistry_wiseman.csv",
                         header = T)
    acdom.data <- read.csv("./data/ag_long_wiseman.csv",
                           header = T)
    anap.data <- read.csv("./data/ad_long_wiseman.csv",
                          header = T)
    
    chl.sample <- bgc.data$Chl
    
    acdom.440.sample <- acdom.data$ag[acdom.data$wavelength == "440"]
    
    anap.440.sample <- anap.data$ad[anap.data$wavelength == "440"]
    
    adg.sample = acdom.440.sample + anap.440.sample
    
    #Truncate the vector of prior parameters
    if (!all(is.na(truncate_chl))) {
      chl.sample = chl.sample[chl.sample >= truncate_chl[1] & chl.sample <= truncate_chl[2]]
      
    }
    
    if (!all(is.na(truncate_adg))) {
      adg.sample = adg.sample[adg.sample >= truncate_adg[1] & adg.sample <= truncate_adg[2]]
      
    }
    
    if (!is.na(sample_count)) {
      set.seed(123)
      chl.sample = sample(chl.sample, size = sample_count)
      adg.sample = sample(adg.sample, size = sample_count)
    }
    
    #try to fit the user-defined distribution
    fit.chl.norm <- fitdistrplus::fitdist(chl.sample, distrib.fit) 
    
    fit.acdm440.norm <- fitdistrplus::fitdist(adg.sample,distrib.fit) 
    
  }
  
  if (use.nomad.prior == TRUE) {
    
    #Load NOMAD in situ Prior data
    bgc.data <- read.csv("./data/nomad_dataset_simplified.csv",
                         header = T, sep = ",")
    
    
    chl.sample <- bgc.data$chl
    chl.sample = chl.sample[!chl.sample %in% -999]
    
    
    acdom.440.sample <- bgc.data$ag443
    acdom.440.sample = acdom.440.sample[!acdom.440.sample %in% -999]
    
    
    anap.440.sample <- bgc.data$ad443
    anap.440.sample = anap.440.sample[!anap.440.sample %in% -999]
    
    adg.sample = acdom.440.sample + anap.440.sample
    
    bbp.sample = bgc.data$bb555
    bbp.sample = bbp.sample[!bbp.sample %in% -999]
    
    #Truncate the vector of prior parameters
    if (all(!is.na(truncate_chl))) {
      chl.sample = chl.sample[chl.sample >= truncate_chl[1] & chl.sample <= truncate_chl[2]]
      
    }
    
    if (all(!is.na(truncate_adg))) {
      adg.sample = adg.sample[adg.sample >= truncate_adg[1] & adg.sample <= truncate_adg[2]]
      
    }
    
    if (all(!is.na(truncate_bbp))) {
      bbp.sample = bbp.sample[bbp.sample >= truncate_bbp[1] & bbp.sample <= truncate_bbp[2]]
      
    }
    
    if (!is.na(sample_count)) {
      set.seed(123)
      chl.sample = sample(chl.sample, size = sample_count)
      adg.sample = sample(adg.sample, size = sample_count)
      bbp.sample = sample(bbp.sample, size = sample_count)
    }
    
    #try to fit the user-defined distribution
    fit.chl.norm <- fitdistrplus::fitdist(chl.sample, distrib.fit) 
    
    fit.acdm440.norm <- fitdistrplus::fitdist(adg.sample,distrib.fit) 
    
    fit.bbp550.norm <- fitdistrplus::fitdist(bbp.sample,distrib.fit)
    
    
  }
  if (plot.diag == TRUE && use.wise.prior == TRUE) {
    
    plot(fit.chl.norm)
    plot(fit.acdm440.norm) 
    
  } else {
    plot(fit.chl.norm)
    plot(fit.acdm440.norm)
    plot(fit.bbp550.norm)
  }
  
  if (plot.diag == TRUE && use.wise.prior == TRUE) {
    
    return(list("fit_chl"=fit.chl.norm, "obs_chl"= chl.sample,
                "fit_acdm440"=fit.acdm440.norm, "obs_acdm440"= adg.sample,
                #"fit_bbp555"=fit.bbp550.norm, "obs_bbp550" = bbp.sample,
                "distribution.fitted" = distrib.fit
    ))
    
  } else {
    return(list("fit_chl"=fit.chl.norm, "obs_chl"= chl.sample,
                "fit_acdm440"=fit.acdm440.norm, "obs_acdm440"= adg.sample,
                "fit_bbp555"=fit.bbp550.norm, "obs_bbp550" = bbp.sample,
                "distribution.fitted" = distrib.fit
    ))
  }
  
  
  
}

#-------------------------------------------------------------------------
#Create Prior density function
#-------------------------------------------------------------------------
prior_bayes = function(param, pop_sd = pop.sd, verbose = F){
  chl = param[1]
  acdom.440 = param[2]
  anap.440 = param[3]
  x.sd = param[4]
  chl.prior = dweibull(x=chl,shape = fit.chl.norm$estimate[1], 
                       scale =fit.chl.norm$estimate[2] , log = T)
  
  acdom440.prior = dweibull(x=acdom.440,shape = fit.acdom440.norm$estimate[1], 
                            scale =fit.acdom440.norm$estimate[2] , log = T)
  
  anap440.prior = dweibull(x=anap.440,shape = fit.anap440.norm$estimate[1], 
                           scale =fit.anap440.norm$estimate[2] , log = T)
  
  lklhood.prior = dunif(x = x.sd, min = 0.000001, max  = 0.001, log = T)
  
  if (verbose == T) {
    if (pop_sd == TRUE) {
      print("Prior of forward model noise is not fitted")
      return(chl.prior+acdom440.prior+anap440.prior)
      
    } else {
      print("Prior of forward model noise is fitted")
      return(chl.prior+acdom440.prior+anap440.prior+lklhood.prior)
    }
  }
  
  
  #
}

#-------------------------------------------------------------------------
#Create Prior sampling function
#-------------------------------------------------------------------------
sampler = function(n=1, pop_sd = pop.sd, verbose=F ){
  
  chl.prior = rweibull(n, shape = fit.chl.norm$estimate[1], 
                       scale =fit.chl.norm$estimate[2])
  
  acdom440.prior = rweibull(n,shape = fit.acdom440.norm$estimate[1], 
                            scale =fit.acdom440.norm$estimate[2])
  
  anap440.prior = rweibull(n,shape = fit.anap440.norm$estimate[1], 
                           scale =fit.anap440.norm$estimate[2])
  
  lklhood.prior = runif(n, min = 0.0001, max = 0.01)
  
  if (verbose == TRUE) {
    if (pop_sd == "known") {
      
      print("Prior of forward model noise is not fitted")
      return(cbind(chl.prior,acdom440.prior,anap440.prior))
      
    } else {
      
      print("Prior of forward model noise is fitted")
      return(cbind(chl.prior,acdom440.prior,anap440.prior,lklhood.prior))
    }
  }
  
}

#-------------------------------------------------------------------------
#Create likelihood function
#-------------------------------------------------------------------------
# initial = par0
# initial_rb_length = length(initial[5:(length(initial)-1)])
NLL_unconstr = function(pars) {
  
  # Values predicted by the forward model
  if (sa.model == "am03") {
    Gpred = Saber_forward_fast(
      use_true_IOPs = F, 
      
      chl = pars[1], 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = pars[4],
      rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
      
      
      Rrs_input_for_slope = obsdata,
      
      slope.parametric = auto_spectral_slope,
      
      
      use_manual_slope =manual_spectral_slope,
      manual_slope =  manual_spectral_slope_vals,
      
      verbose = F, wavelength = wavelength
    )
    
  } else {
    Gpred = lee_forward_fast(
      use_true_IOPs = F, 
      
      chl = pars[1], 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = pars[4],
      rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
      
      
      Rrs_input_for_slope = obsdata,
      
      slope.parametric = auto_spectral_slope,
      
      
      use_manual_slope =manual_spectral_slope,
      manual_slope =  manual_spectral_slope_vals,
      
      verbose = F, wavelength = wavelength
    )
  }
  
  # Negative log-likelihood
  smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                     log = TRUE), na.rm = T)
  return(smull)
}


#-------------------------------------------------------------------------
#Create a single Bayes Function for the MCMC sampling
#-------------------------------------------------------------------------

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

runBayes <- function(obsdata, rrs_type, max_par, min_par, 
                     constrain_config = c("const_bbp" = F, "cost_bgc"=F, "cost_IOP"=F),
                     qaa_slope, manual_slope, 
                     manual_slope_val= c("s_g"=0.014, "s_d"=0.003, "gamma"=0.5),
                     iter_count, sampler_mcmc,
                      wavelngth_sim, sa_model, hybrid_mode){
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### EMBRACE THE RANDOMNESS #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  cat(paste0("\033[0;34m MODEL CONSTRAINTS ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Backscatter constrain: [chl]=",constrain_config[1],"\033[0m","\n"))
  cat(paste0("\033[0;32m","Biogeochemical variable constrain: [chl]=",constrain_config[2],"\033[0m","\n"))
  cat(paste0("\033[0;32m","IOP constrain: [chl]=",constrain_config[3],"\033[0m","\n"))
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m SPECTRAL SLOPE SETTINGS ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope for QAA : [chl]=",qaa_slope,"\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope from user: [chl]=",manual_slope,"\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope values: [chl]=",manual_slope_val,"\033[0m","\n"))
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m ANCILIARY ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Rrs type: ",rrs_type,"\033[0m","\n"))
  if (sa_model == "am03") {
    cat(paste0("\033[0;32m","Paramterization for forward model: Albert & Mobley (2003) \033[0m","\n"))
    
  } else {
    cat(paste0("\033[0;32m","Paramterization for forward model: Lee et al., (1998) \033[0m","\n"))
    
  }
  cat(paste0("\033[0;32m","MCMC iteration count:",iter_count,"\033[0m","\n"))
  cat(paste0("\033[0;32m","MCMC Sampler:",sampler_mcmc,"\033[0m","\n"))
  cat(paste0("\033[0;32m","HYBRID mode: ",hybrid_mode,"\033[0m","\n"))
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  Sys.sleep(1)
  
  LL_unconstr_bayes = function(pars) {
    
    # Values predicted by the forward model
    if (sa_model == "am03") {
      Gpred = Saber_forward_fast(
        use_true_IOPs = F, 
        
        chl = pars[1], 
        a_dg = pars[2],
        bbp.550 = pars[3],
        
        z = pars[4],
        rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
        
        
        Rrs_input_for_slope = obsdata,
        
        slope.parametric = qaa_slope,
        
        
        use_manual_slope =manual_slope,
        manual_slope =  manual_spectral_slope_vals,
        
        verbose = F, wavelength = wavelngth_sim
      )
      
    } else {
      Gpred = lee_forward_fast(
        use_true_IOPs = F, 
        
        chl = pars[1], 
        a_dg = pars[2],
        bbp.550 = pars[3],
        
        z = pars[4],
        rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
        
        
        Rrs_input_for_slope = obsdata,
        
        slope.parametric = qaa_slope,
        
        
        use_manual_slope =manual_slope,
        manual_slope =  manual_spectral_slope_vals,
        
        verbose = F, wavelength = wavelngth_sim
      )
    }
    
    # Negative log-likelihood
    smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                      log = TRUE), na.rm = T)
    return(smull)
  }
  
  #Instantiate inversion objective function and the optimization scheme 
  
  inv_bound = create_init_bound(rrs_inv_type = rrs_type, manual_par0 = T, 
                                constrain.bbp = model_config[1], 
                                constrain.shallow.bgc = model_config[2], 
                                constrain.shallow.iop = model_config[3], 
                                pop.sd =  "unknown",
                                
                                init_par = rowMeans(matrix(c(max_par,min_par), 
                                                           ncol=2)),
                                upper_par = max_par,
                                lower_par = min_par
  )
  
  par0_bayes = inv_bound$par0 ; upper.bound_bayes = inv_bound$upper.bound  
  lower.bound_bayes= inv_bound$lower.bound
  
  cat(paste0("\033[0;32m","MINIMA values are: [chl]=",lower.bound_bayes[1], 
             ",adg(440)=",lower.bound_bayes[2], 
             ",bbp(555)=", lower.bound_bayes[3], ", zB=", 
             lower.bound_bayes[4],",RB={", 
             toString(as.numeric(lower.bound_bayes[5:(length(lower.bound_bayes)-1)])),
             "},population_sigma=", lower.bound_bayes[length(lower.bound_bayes)],"\033[0m","\n"))
  
  cat(paste0("\033[0;32m","MAXIMA values are: [chl]=",upper.bound_bayes[1], 
             ",adg(440)=",upper.bound_bayes[2], 
             ",bbp(555)=", upper.bound_bayes[3], ", zB=", 
             upper.bound_bayes[4],",RB={", 
             toString(as.numeric(upper.bound_bayes[5:(length(upper.bound_bayes)-1)])),
             "},population_sigma=", upper.bound_bayes[length(upper.bound_bayes)],"\033[0m","\n"))
  
  
  if (hybrid_mode == TRUE) {
    
    cat(paste0("\033[0;34m ENTERING HYBRID MODE ::::\033[0m","\n"))
    
    inverse_output <- quiet(suppressWarnings(solve.objective.inverse.shallow.final.fast(
      
      wave = wavelngth_sim, 
      
      constrain.shallow.iop = model_config[3],  
      
      unconstrained = T,
      
      constrained.bgc = model_config[2], 
      constrain.bgc.value = NA,
      
      
      initial = as.numeric(par0_bayes), 
      obsdata = as.numeric(obsdata),
      
      auto_spectral_slope = qaa_slope,
      manual_spectral_slope = manual_slope, 
      
      manual_spectral_slope_vals = manual_slope_val,
      
      sa.model = sa_model, 
      
      obj.fn = "log-LL",
      method.opt = methods.opt[4],
      
      lower.b = lower.bound_bayes,
      upper.b = upper.bound_bayes, 
    )))
    
    # Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
    #                                inverse_output[[1]]$`sd(+/-)`)
    
    bayessetup <- createBayesianSetup(prior = NULL,
                                      likelihood = LL_unconstr_bayes,
                                      lower = lower.bound,#[c(1:3,length(lower.bound))] 
                                        best = inverse_output[[1]]$estimates, 
                                      upper = upper.bound#[c(1:3,length(lower.bound))]
                                      ,
                                      names = c("chl","adg443","bbp555", 
                                                "H",  rep(paste0("fa", (1:rb_count))), 
                                                "pop_sd"),
                                      parallel = F)
    
  } else {
    
    bayessetup <- createBayesianSetup(prior = NULL,
                                      likelihood = LL_unconstr_bayes,
                                      lower = lower.bound#[c(1:3,length(lower.bound))]
                                      , 
                                      upper = upper.bound#[c(1:3,length(lower.bound))]
                                      ,
                                      names = c("chl","adg443","bbp555", 
                                                "H",  rep(paste0("fa", (1:rb_count))), 
                                                "pop_sd"),
                                      parallel = F)
    
  }
  
  cat(paste0("\033[0;34m @@@@@@@@@@ TEST Bayesian Set-UP @@@@@@@@@ \033[0m","\n"))
  checkBayesianSetup(bayessetup) #Test if the setup is initiated for theta pars
  
  max_iter = iter_count ; burn_in = max_iter/4; chain_burn = burn_in/3
  
  settings = list(iterations = max_iter, message = TRUE, nrChains = 1 #Set MCMC config
                  , burnin=burn_in
  )
  
  #Run MCMC
  cat(paste0("\033[0;34m @@@@@@@@@@ MCMC START @@@@@@@@@ \033[0m","\n"))
  out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = sampler_mcmc )
  
  sd_c1 = apply(out[["chain"]][[1]][-1,1:8], 2, sd)
  sd_c2 = apply(out[["chain"]][[2]][-1,1:8], 2, sd)
  sd_c3 = apply(out[["chain"]][[3]][-1,1:8], 2, sd)
  
  sd_bayes = rowMeans(matrix(c(sd_c1,sd_c2,sd_c3), 
                             ncol=3))
  MAP_bayes = MAP(out, start = burn_in
                  )
  bayes_output = c(MAP_bayes$parametersMAP, sd_bayes)
  cat(paste0("\033[0;34m @@@@@@@@@@ MCMC FINISHED @@@@@@@@@ \033[0m","\n"))
  summary(out, start = burn_in)
  
  return(bayes_output)
  
}


runBayes(obsdata = obsdata, rrs_type = "shallow", max_par = c(50, 5, 0.01, 10, 1,1,1, 10), 
         min_par = c(0.5, 0.1, 0.001, 0.5, 0,0,0, 0.0001), constrain_config = c(F, F, F), 
         qaa_slope = F, manual_slope = F, iter_count = 20000, sampler_mcmc = samplerlist[6], 
         wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F )

param_vec[test_idx,]











