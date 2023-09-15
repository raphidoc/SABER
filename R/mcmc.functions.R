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

#Create Prior density function
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

#Create Prior sampling function
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


#Create likelihood function

# ll <-function(param){
#   Gpred = Saber_forward(chl = param[1], acdom440 = param[2],
#                         anap440 =param[3],
#                         #bbp.550 = HL.deep.iop$bbp550[j],
#                         bbp.550 = Fit.input$bbp.550,
#                         verbose=F ,realdata = obsdata)
#   
#   # smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.0006327431, log = TRUE),
#   #             na.rm = T)
#   # 
#   # smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.011142, log = TRUE),
#   #             na.rm = T)
#   
#   smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = param[4], 
#                     log = TRUE),na.rm = T) #10000 is the scaling factor to avoid calculation  
#                                            #of very small numbers
#   
#   #smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE))
#   
#   return(smull)
# }


initial = par0
initial_rb_length = length(initial[5:(length(initial)-1)])
NLL_unconstr = function(pars) {
  
  # Values predicted by the forward model for single RUN
  if (sa.model == "am03") {
    Gpred = Saber_forward_fast(
      use_true_IOPs = F, 
      #a_non_water_path = IOP_files[idx_a],
      #bb_non_water_path = IOP_files[idx_bb],
      
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
    Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = bbp.550,
                        z = pars[4], rb.fraction = pars[5:(length(pars)-1)],
                        verbose = F, realdata = data, plot = F)
  }
  
  # Negative log-likelihood
  smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                     log = TRUE), na.rm = T)
  return(smull)
}
