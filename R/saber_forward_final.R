
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
# Function to read the configuration file for forward modeling
#=======================================================================
read_fwd_param <- function(fwd_config_file){
  
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