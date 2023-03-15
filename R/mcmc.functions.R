create.prior.data <- function(use.ioccg.prior=F, use.wise.prior=F, use.nomad.prior=T,
                              plot.diag=F){
  
  if (use.ioccg.prior == TRUE) {
    
    Hl.deep.iop = read_IOCCG_data()
    
    fit.chl.norm <- fitdistrplus::fitdist(HL.deep.iop$chl, "weibull") #try to fit dist
    
    
    fit.acdom440.norm <- fitdistrplus::fitdist(HL.deep.iop$acdom440, "weibull") #try to fit dist
    
    
    fit.anap440.norm <- fitdistrplus::fitdist(HL.deep.iop$anap440, "weibull") #try to fit dist
    
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
    #hist(acdom.440.sample, probability = T)
    
    anap.440.sample <- anap.data$ad[anap.data$wavelength == "440"]
    #hist(anap.440.sample, probability = T)
    
    fit.chl.norm <- fitdistrplus::fitdist(chl.sample, "weibull") #try to fit 
    #plot(fit.chl.norm) #Infer results of the fit
    
    fit.acdom440.norm <- fitdistrplus::fitdist(acdom.440.sample, "weibull") #try to fit 
    
    
    
    fit.anap440.norm <- fitdistrplus::fitdist(anap.440.sample, "weibull") #try to fit 
    #plot(fit.anap440.norm) #Infer results of the fit
    
    
  }
  
  if (use.nomad.prior == TRUE) {
    
    #Load NOMAD in situ Prior data
    bgc.data <- read.csv("./data/nomad_dataset_simplified.csv",
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
    
    
  }
  if (plot.diag == TRUE) {
    
    plot(fit.acdom440.norm) #Infer results of the fit
    plot(fit.chl.norm) #Infer results of the fit
    plot(fit.anap440.norm) #Infer results of the fit
    
  }
  
  return(list("fit.chl"=fit.chl.norm,
  "fit.acdom440"=fit.acdom440.norm,
  "fit.anap440"=fit.anap440.norm))
  
}

#Create Prior density function
prior = function(param, pop_sd = pop.sd, verbose = F){
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
ll <-function(param){
  Gpred = Saber_forward(chl = param[1], acdom440 = param[2],
                        anap440 =param[3],
                        #bbp.550 = HL.deep.iop$bbp550[j],
                        bbp.550 = Fit.input$bbp.550,
                        verbose=F ,realdata = obsdata)
  
  # smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.0006327431, log = TRUE),
  #             na.rm = T)
  # 
  # smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.011142, log = TRUE),
  #             na.rm = T)
  
  smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = param[4], 
                    log = TRUE),na.rm = T) #10000 is the scaling factor to avoid calculation  
                                           #of very small numbers
  
  #smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE))
  
  return(smull)
}
