
source("./R/saber_forward_parametric_conc_wise.R") #FINAL VERSION
source("./R/SABER_forward_fast.R") #FASTEST VERSION

#Inelastic Scattering Models
source("./R/SABER.sicf.R")
source("./R/SABER.fdom.R")
source("./R/inelastic.wavelength.redistribution.R")

#QAAv5 for IOP and Spectral slopes
source("./R/QAA_v5.R")

#Inverse Am03 model with residual function
#source("./R/saber_inverse.R")

#Sun Earth Geometry and Downwelling Irradiance 
source("./R/snell_law.R")
source("./R/gregg_carder_ed.R")

#Spectral Bottom Reflectance Retrieval
source("./R/retrive.rb.spectral.wise.R")
source("./R/generate.rb.spectral.R")
source("./R/read.surface.IOPs.wise.R")
source("./R/retrieve.rb.spectral.minimal.R")

#Gradient Based Inverse Functions
source("./R/solve.objective.inverse.final.R") 
source("./R/solve.objective.inverse_fast.R") #FINAL VERSION

#MCMC Based Inverse Functions
source("./R/mcmc.functions.R")
source("./R/saber.inversion.vector.parallel.R")
source("./R/mcmc_bayes_parallel.R")
#-------------------------------------------------------------------------------------------
library(dplyr)
library(readxl)
library(stats4)
library(MASS)
library(dglm)
library(fitdistrplus)
library(Riops)
library(Cops)
library(ggplot2)
library(rho)
library(marqLevAlg)
library(BayesianTools)
library(coda)
library(stats)
library(Deriv)
library(mcmc)
library(plot3D)
library(limSolve)
library(nloptr)
library(minpack.lm)
library(pracma)
library(alabama)
library(data.table)
library(ggalt)
library(ggExtra)
library(suncalc)
library(viridis)

#=======================================================================
# Function to read the configuration file for inversion
#=======================================================================
read_inv_param <- function(inv_config_file){
  
  inv_config = read.table(file = inv_config_file, sep=",", header = T,colClasses = "character")
  rrs_type = as.character(trimws(inv_config$param_vals[1]))
  
  mode = as.character(trimws(inv_config$param_vals[2]))
  
  param_name = as.character(trimws((unlist((strsplit(inv_config$param_vals[3],";"))))))
  min_par = as.numeric((unlist((strsplit(inv_config$param_vals[4],";")))))
  max_par = as.numeric((unlist((strsplit(inv_config$param_vals[5],";")))))
  init_par = as.numeric((unlist((strsplit(inv_config$param_vals[6],";")))))
  
  constrain_config = as.logical(trimws((unlist((strsplit(inv_config$param_vals[7],";"))))))
  names(constrain_config) = c("const_bbp", "const_bgc", "const_IOP")
  constrain_bbp_val = as.numeric(inv_config$param_vals[8])
  constrain_bgc_val = as.numeric(trimws((unlist((strsplit(inv_config$param_vals[9],";"))))))
  names(constrain_bgc_val) = c("chl", "adg443", "bbp555")
  constrain_iop = as.character(trimws((unlist((strsplit(inv_config$param_vals[10],";"))))))
  
  qaa_mode = as.logical(trimws(inv_config$param_vals[11]))
  qaa_prefit = as.logical(trimws(inv_config$param_vals[12]))
  
  auto_slope = as.logical(trimws(inv_config$param_vals[13]))
  manual_slope  =  as.logical(trimws(inv_config$param_vals[14]))
  manual_slope_vals = as.numeric((unlist((strsplit(inv_config$param_vals[15],";")))))
  names(manual_slope_vals) = c("s_g", "s_d", "gamma")
  
  sa_model = as.character(trimws((unlist((strsplit(inv_config$param_vals[16],";"))))))
  obj_fn = as.character(trimws((unlist((strsplit(inv_config$param_vals[17],";"))))))
  opt_method = as.character(trimws((unlist((strsplit(inv_config$param_vals[18],";"))))))
  
  iter_count = as.numeric(trimws(inv_config$param_vals[19]))
  sampler_mcmc = as.character(trimws(inv_config$param_vals[20]))
  hybrid_mode = as.logical(trimws(inv_config$param_vals[21]))
  
  config_list = list(rrs_type, mode, param_name, min_par, max_par, init_par, constrain_config,
                     constrain_bbp_val, constrain_bgc_val, constrain_iop, qaa_mode,
                     qaa_prefit, auto_slope, manual_slope, manual_slope_vals, 
                     sa_model,obj_fn, opt_method, iter_count, sampler_mcmc, hybrid_mode )
  names(config_list) = inv_config$param_config
  return(config_list)

}

#test_param = read_inv_param(inv_config_file = "./data/inversion_config.txt")

#=======================================================================
# Function to create initial values and bounds for inversion
#=======================================================================
rb_count = length(fA.set)
create_init_bound <- function(pop.sd = "unknown", 
                              constrain.bbp= F,
                              constrain.shallow.bgc = F , constrain.shallow.iop = F,
                              manual_par0 = T, 
                              
                              init_par = c(chl = 2, adg440 = 0.8,
                                           bbp = 0.005,z = 2.5, 
                                           rb_c = rep(0.5, rb_count),
                                           population.sd = 0.1),
                              
                              upper_par = c(chl = 30, adg440 = 5,
                                            bbp = 0.01,z = 10, 
                                            rb_c = rep(1, rb_count),
                                            population.sd = 10),
                              
                              lower_par = c(chl = 1, adg440 = 0.2,
                                            bbp = 0.001,z = 0.5, 
                                            rb_c = rep(0, rb_count),
                                            population.sd = 1e-05),
                              
                              rrs_inv_type = "deep"){
  
  #-------------------------------------
  # Initial values for deep water 
  #-------------------------------------
  if (pop.sd == "unknown" & rrs_inv_type == "deep" & manual_par0 == TRUE & constrain.bbp == "TRUE") {# <<MANUAL INPUT>>
    par0 = c(chl = init_par[1],
             adg = init_par[2], 
             #bbp550 =0.005,
             population.sd = init_par[length(init_par)]) #@@  sigma UNKNOWN
  }
  
  if (pop.sd == "unknown" & rrs_inv_type == "deep" & manual_par0 == TRUE & constrain.bbp == "FALSE") {# <<MANUAL INPUT>>
    par0 = c(chl = init_par[1],
             adg = init_par[2], 
             bbp550 =init_par[3],
             population.sd = init_par[length(init_par)]) #@@  sigma UNKNOWN
  }
  
  #Set bounds (Instead of flat multiplier, implement range as a function of parameter sensitivity)
  if (rrs_inv_type == "deep" & constrain.bbp == TRUE & manual_par0 == FALSE) {
    upper.bound <- c(5*init_par[1:2], upper_par[length(upper_par)])  
    lower.bound <- c(0.2*init_par[1:2], lower_par[length(lower_par)])
    
  } else {
    
    if (rrs_inv_type == "deep" & manual_par0 == T) {
      if (constrain.bbp == FALSE) {
        upper.bound <- upper_par[c(1:3,length(upper_par))]  
        lower.bound <- lower_par[c(1:3,length(lower_par))]
      } else {
        upper.bound <- upper_par[c(1:2,length(upper_par))]  
        lower.bound <- lower_par[c(1:2,length(lower_par))]
      }
      
    }
    
  }  
  
  #--------------------------------------
  # Initial values for shallow water
  #--------------------------------------
  rb_count = length(fA.set)
  rep(0.5, rb_count)
  rb_c = paste0("rb", (1:rb_count))
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.bgc == "TRUE") { # <<MANUAL INPUT>>
    
    par0 = init_par[-(1:3)]
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par[-(1:3)]
      upper.bound = upper_par[-(1:3)]
      
    } else {
      
      lower.bound <- c(z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c(z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
    
  }
  
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.iop == "TRUE") {
    par0 = c(z = init_par[4], 
             rb_c = init_par[5:7],
             population.sd = init_par[length(init_par)])
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par[4:length(lower_par)]
      upper.bound = upper_par[4:length(upper_par)]
      
    } else {
      
      lower.bound <- c(z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c(z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
    
  }
  
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.bgc == "FALSE" & 
      constrain.shallow.iop == "FALSE") { 
    
    par0 = init_par
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par
      upper.bound = upper_par
      
    } else {
      
      #Autoscale Intital values from pre-Ft
      increament.scale <- 1
      
      lower.bound <- c((par0[1:3] - 0.8*par0[1:3]),z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c((par0[1:3] + 5*par0[1:3]),z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
  }
  
  return(list("par0"=par0, "upper.bound"=upper.bound, "lower.bound"=lower.bound))
}

# inv_bound = create_init_bound(rrs_inv_type = rrs_type, manual_par0 = T, 
#                               constrain.bbp = F, 
#                               constrain.shallow.bgc = T, 
#                               constrain.shallow.iop = F, 
#                               pop.sd =  "unknown",
#                               
#                               init_par = rowMeans(matrix(c(max_par,min_par), 
#                                                          ncol=2)),
#                               upper_par = max_par,
#                               lower_par = min_par
# )
# 
# par0 = inv_bound$par0 ; upper_bound = inv_bound$upper.bound  
# lower_bound= inv_bound$lower.bound

#=======================================================================
# A single Function to perform the inversion as per the configuration
#=======================================================================
run_inverse_mode <- function(rrs_input, 
                             wavelength_input, 
                             inv_config_file = "./data/inversion_config.txt"){
  
  config_list = read_inv_param(inv_config_file)
  
  print("READING CONFIGURATION FILE...")
  
  if (config_list$` rrs_type` == "deep") {
    type_Rrs_below = "deep"
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    if (config_list$` mode` == "use_grad") {
    inv_output = inverse_runGrad(rrs_type = config_list$` rrs_type`, obsdata = rrs_input, 
                                 max_par = config_list$` max_par`, 
                                 min_par = config_list$` min_par`, 
                                 init_par = config_list$` init_par`, 
                                 param_lab = config_list$` param_name`, 
                                 constrain_config = config_list$` constrain_conf`, 
                                 bgc_const_val = config_list$` constrain_bgc_val`, 
                                 bbp_const_val = config_list$` constrain_bbp_val`,
                                 iop_const_path =  config_list$` constrain_iop`,
                                 qaa_prefit = config_list$` qaa_prefit`, 
                                 QAA_mode = config_list$` qaa_mode`, 
                                 qaa_slope = config_list$` auto_slope`, 
                                 manual_slope = config_list$` manual_slope`, 
                                 manual_slope_vals = config_list$` manual_slope_val`, 
                                 wavelength_sim = wavelength_input, 
                                 sa_model = config_list$` sa_model`, 
                                 obj_fn = config_list$` obj_fn`, 
                                 opt_method = config_list$` opt_method`)
    }
    
    if (config_list$` mode` == "use_mcmc") {
      inv_output = inverse_runBayes(obsdata = rrs_input, 
                                    rrs_type =  type_Rrs_below, 
                                    max_par = config_list$` max_par`, 
                                    min_par = config_list$` min_par`, 
                                    param_lab = config_list$` param_name`, 
                                    constrain_config = config_list$` constrain_conf`, 
                                    bbp_const_val = config_list$` constrain_bbp_val`, 
                                    bgc_const_val = config_list$` constrain_bgc_val`, 
                                    iop_const_path = config_list$` constrain_iop`,
                                    qaa_slope = config_list$` auto_slope`, 
                                    manual_slope = config_list$` manual_slope`, 
                                    manual_slope_val = config_list$` manual_slope_val`, 
                                    iter_count = config_list$` iter_count`, 
                                    sampler_mcmc = config_list$` sampler_mcmc`, 
                                    wavelngth_sim = wavelength_input, 
                                    sa_model = config_list$` sa_model`, 
                                    hybrid_mode = config_list$` hybrid_mode`, plot_rrs = T)
    }
    
    
  }
  
  if (config_list$` rrs_type` == "shallow") {
    type_Rrs_below = "shallow"
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    if (config_list$` mode` == "use_grad") {
        
        inv_output = inverse_runGrad(rrs_type = config_list$` rrs_type`, 
                                     obsdata = rrs_input, max_par = config_list$` max_par`, 
                                     min_par = config_list$` min_par`, 
                                     init_par = config_list$` init_par`, 
                                     param_lab = config_list$` param_name`, 
                                     constrain_config = config_list$` constrain_conf`, 
                                     bgc_const_val = config_list$` constrain_bgc_val`,
                                     bbp_const_val = config_list$` constrain_bbp_val`,
                                     iop_const_path =  config_list$` constrain_iop`,
                                     qaa_prefit = config_list$` qaa_prefit`, 
                                     QAA_mode = config_list$` qaa_mode`, 
                                     qaa_slope = config_list$` auto_slope`, 
                                     manual_slope = config_list$` manual_slope`, 
                                     manual_slope_vals = config_list$` manual_slope_val`, 
                                     wavelength_sim = wavelength_input, 
                                     sa_model = config_list$` sa_model`, 
                                     obj_fn = config_list$` obj_fn`, 
                                     opt_method = config_list$` opt_method`)
      
    }
    
    if (config_list$` mode` == "use_mcmc") {
      type_Rrs_below = config_list$` rrs_type`
      assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
      
      inv_output = inverse_runBayes(obsdata = rrs_input, 
                                    rrs_type =  type_Rrs_below, 
                                    max_par = config_list$` max_par`, 
                                    min_par = config_list$` min_par`, 
                                    param_lab = config_list$` param_name`, 
                                    constrain_config = config_list$` constrain_conf`,
                                    bbp_const_val = config_list$` constrain_bbp_val`, 
                                    bgc_const_val = config_list$` constrain_bgc_val`, 
                                    iop_const_path = config_list$` constrain_iop`,
                                    qaa_slope = config_list$` auto_slope`, 
                                    manual_slope = config_list$` manual_slope`, 
                                    manual_slope_val = config_list$` manual_slope_val`, 
                                    iter_count = config_list$` iter_count`, 
                                    sampler_mcmc = config_list$` sampler_mcmc`, 
                                    wavelngth_sim = wavelength_input, 
                                    sa_model = config_list$` sa_model`, 
                                    hybrid_mode = config_list$` hybrid_mode`, plot_rrs = T)
      
    }
    
  }
  
  
  return(inv_output)
}

#=======================================================================
# Test the inversion function
#=======================================================================
run_inverse_mode(rrs_input = obsdata #+ rnorm.trunc(n = length(obsdata), mean = 0, 
                                             #sd = sd(obsdata), min = 0)*0.05 
                      , wavelength_input = wavelength)

run_inverse_mode(rrs_input = as.numeric(rrs.forward.SABER[1,]), 
                 wavelength_input = wavelength)
param_vec[100,]


