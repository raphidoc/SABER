#-----------------------------------------------------------------------------------------
#Vectorize the inversion for multiple Rrs observation
#-----------------------------------------------------------------------------------------

#===================================================================
### Deep water unconstrained inversion
#===================================================================

# Unconstrained inversion for deep water (with QAA support)
doOptimization_deep_unconst <- function(obsdata, init_par, max_par, min_par,
                                        wavelength_sim, 
                                        qaa_prefit, QAA_mode, 
                                        qaa_slope, manual_slope, manual_slope_vals,
                                        sa_model, obj_fn, opt_method) {
  
  #Assign the inversion vector of rrs
  rrs_inverse_input = obsdata
  
  if (QAA_mode == F) {
    
    if (qaa_prefit == TRUE) {
      print("QAA PREFIT : ON")
      qaa_output = QAA.v5(waves = wavelength_sim, Rrs = rrs_inverse_input)
      init_par = c(chl = qaa_output$chl, adg440 =qaa_output$a_dg_443, 
               bbp550 = qaa_output$b_bp_555, "pop_sd" = init_par[length(init_par)])
      
      if (any(init_par == Inf) | any(is.na(init_par))) {
        print("Prefit failed, default values wil be used")
        init_par = init_par
        
        max_par <- max_par  
        min_par <- min_par
      } else {
        print(paste0("prefit values are computed as chl: ", signif(init_par[1], digits = 3), 
                     " adg443: ", signif(init_par[2], digits = 3), 
                     " bbp555: ",signif(init_par[3], digits = 3), 
                     " pop_sd: ",signif(init_par[4], digits = 3)))
        
        max_par <- c(init_par[1:2] + 0.15*init_par[1:2], max_par[3:4])  
        min_par <- c(init_par[1:2] - 0.15*init_par[1:2], min_par[3:4])
      }
      
    } else {
      print("QAA PREFIT : OFF")
      init_par = init_par
      print(paste0("prefit values are computed as chl: ", signif(init_par[1], digits = 3), 
                   " adg443: ", signif(init_par[2], digits = 3), 
                   " bbp555: ",signif(init_par[3], digits = 3), 
                   " pop_sd: ",signif(init_par[4], digits = 3)))
      max_par <- max_par  
      min_par <- min_par
    }
    
    inverse_output_deep <- suppressWarnings(solve.objective.inverse.deep.final.fast(
      
      
      wave = wavelength_sim,
      initial = init_par, 
      bbp.constrain  = F,
      bbp.constrain.value =  NA,
      
      obsdata = as.numeric(rrs_inverse_input),
      
      sa.model = sa_model, 
      obj.fn = obj_fn,
      
      auto_spectral_slope = qaa_slope,
      manual_spectral_slope = manual_slope, 
      
      manual_spectral_slope_vals = manual_slope_vals,
      
      method.opt = opt_method,
      lower.b = min_par,
      upper.b = max_par 
      #,batch = FALSE, pop.sd = FALSE
    ))
    
    Fit.optimized.ssobj.batch <- c(inverse_output_deep[[1]]$estimates,
                                   inverse_output_deep[[1]]$`sd(+/-)`)
    
  } else {
    
    
    qaa_output = QAA.v5(waves = wavelength_sim, Rrs = rrs_inverse_input)
    Fit.optimized.ssobj.batch = c(chl = qaa_output$chl, adg443 =qaa_output$a_dg_443, 
                                  bbp555 = qaa_output$b_bp_555, rep(NA,5))
    
    if (any(Fit.optimized.ssobj.batch[1:3] == Inf) | any(is.na(Fit.optimized.ssobj.batch[1:3]))) {
      print("QAA failed, NA values wil be used")
      Fit.optimized.ssobj.batch = rep(NA,8)
    }
    
    return(Fit.optimized.ssobj.batch)
  }
}


#===================================================================
### Deep water constrained inversion
#===================================================================

# Unconstrained inversion for deep water (with QAA support)
doOptimization_deep_const <- function(obsdata, init_par, max_par, min_par,
                                        wavelength_sim, bbp_const_val,
                                        qaa_prefit, QAA_mode, 
                                        qaa_slope, manual_slope, manual_slope_vals,
                                        sa_model, obj_fn, opt_method) {
  
  #Assign the inversion vector of rrs
  rrs_inverse_input = obsdata
  
  if (QAA_mode == F) {
    
    if (qaa_prefit == TRUE) {
      print("QAA PREFIT : ON")
      qaa_output = QAA.v5(waves = wavelength_sim, Rrs = rrs_inverse_input)
      #browser()
      init_par = c(chl = qaa_output$chl, adg440 =qaa_output$a_dg_443, 
                   "pop_sd" = init_par[length(init_par)])
      
      if (any(init_par == Inf) | any(is.na(init_par))) {
        print("Prefit failed, default values wil be used")
        init_par = init_par
        
        max_par <- max_par  
        min_par <- min_par
      } else {
        print(paste0("prefit values are computed as chl: ", signif(init_par[1], digits = 3), 
                     " adg443: ", signif(init_par[2], digits = 3), 
                     #" bbp555: ",signif(init_par[3], digits = 3), 
                     " pop_sd: ",signif(init_par[length(init_par)], digits = 3)))
        
        max_par <- c(init_par[1:2] + 0.15*init_par[1:2], max_par[length(init_par)])  
        min_par <- c(init_par[1:2] - 0.15*init_par[1:2], min_par[length(init_par)])
      }
      
    } else {
      print("QAA PREFIT : OFF")
      init_par = init_par
      print(paste0("prefit values are computed as chl: ", signif(init_par[1], digits = 3), 
                   " adg443: ", signif(init_par[2], digits = 3), 
                   #" bbp555: ",signif(init_par[3], digits = 3), 
                   " pop_sd: ",signif(init_par[4], digits = 3)))
      max_par <- max_par  
      min_par <- min_par
    }
    
    inverse_output_deep <- suppressWarnings(solve.objective.inverse.deep.final.fast(
      
      
      wave = wavelength_sim,
      initial = init_par, 
      bbp.constrain  = T,
      bbp.constrain.value =  bbp_const_val,
      
      obsdata = as.numeric(rrs_inverse_input),
      
      sa.model = sa_model, 
      obj.fn = obj_fn,
      
      auto_spectral_slope = qaa_slope,
      manual_spectral_slope = manual_slope, 
      
      manual_spectral_slope_vals = manual_slope_vals,
      
      method.opt = opt_method,
      lower.b = min_par,
      upper.b = max_par 
      #,batch = FALSE, pop.sd = FALSE
    ))
    
    Fit.optimized.ssobj.batch <- c(inverse_output_deep[[1]]$estimates,
                                   inverse_output_deep[[1]]$`sd(+/-)`)
    
  } else {
    
    
    qaa_output = QAA.v5(waves = wavelength_sim, Rrs = rrs_inverse_input)
    Fit.optimized.ssobj.batch = c(chl = qaa_output$chl, adg443 =qaa_output$a_dg_443, 
                                  rep(NA,2))
    
    if (any(Fit.optimized.ssobj.batch[1:2] == Inf) | any(is.na(Fit.optimized.ssobj.batch[1:2]))) {
      print("QAA failed, NA values wil be used")
      Fit.optimized.ssobj.batch = rep(NA,4)
    }
    
    return(Fit.optimized.ssobj.batch)
  }
}


#===================================================================
# BGC constrained inversion for shallow water
#===================================================================
doOptimization_shallow_IOP_const <- function(obsdata, init_par, max_par, min_par,
                                             bgc_const_val, iop_const_path,
                                             wavelength_sim, 
                                             qaa_prefit, QAA_mode, 
                                             qaa_slope, manual_slope, manual_slope_vals,
                                             sa_model, obj_fn, opt_method) {
  
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
    
    wave = wavelength_sim, 
      
    unconstrained = F, 
    constrain.shallow.iop = T, abs_path = iop_const_path[1], bb_path = iop_const_path[2],
    constrained.bgc  = F, 
    
    constrain.bgc.value = bgc_const_val,
    
    initial = as.numeric(init_par), 
    obsdata = obsdata,
    
    auto_spectral_slope = qaa_slope,
    
    manual_spectral_slope = manual_slope, 
    manual_spectral_slope_vals = manual_slope_vals,
    
    sa.model = sa_model, 
    
    obj.fn = obj_fn,
    method.opt = opt_method,
    
    lower.b = min_par,
    upper.b = max_par
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}


#===================================================================
# BGC constrained inversion for shallow water
#===================================================================
doOptimization_shallow_BGC_const <- function(obsdata, init_par, max_par, min_par,
                                             bgc_const_val,
                                             wavelength_sim, 
                                             qaa_prefit, QAA_mode, 
                                             qaa_slope, manual_slope, manual_slope_vals,
                                             sa_model, obj_fn, opt_method) {
  
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
    
    wave = wavelength_sim,
    
    unconstrained = F, 
    constrain.shallow.iop = F, abs_path = NA, bb_path = NA,
    constrained.bgc  = T, 
    
    constrain.bgc.value = bgc_const_val,
    
    initial = as.numeric(init_par), 
    obsdata = obsdata,
    
    auto_spectral_slope = qaa_slope,
    
    manual_spectral_slope = manual_slope, 
    manual_spectral_slope_vals = manual_slope_vals,
    
    sa.model = sa_model, 
    
    obj.fn = obj_fn,
    method.opt = opt_method,
    
    lower.b = min_par,
    upper.b = max_par 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}




#===================================================================
### Unconstrained inversion for shallow water
#===================================================================
doOptimization_shallow_unconst <- function(obsdata, init_par, max_par, min_par,
                                           bgc_const_val,
                                           wavelength_sim, 
                                           qaa_prefit, QAA_mode, 
                                           qaa_slope, manual_slope, manual_slope_vals,
                                           sa_model, obj_fn, opt_method) {
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
    
    wave = wavelength_sim,
    
    constrain.shallow.iop = F,  
    
    unconstrained = T,
    
    constrained.bgc = F, 
    constrain.bgc.value = bgc_const_val,
    
    
    initial = as.numeric(init_par), 
    obsdata = obsdata,
    
    auto_spectral_slope = qaa_slope,
    
    manual_spectral_slope = manual_slope, 
    manual_spectral_slope_vals = manual_slope_vals,
    
    sa.model = sa_model, 
    
    obj.fn = obj_fn,
    method.opt = opt_method,
    
    lower.b = min_par,
    upper.b = max_par 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}



############################################################################################
# #Instantiate inversion objective function and the optimization scheme 
# pop.sd = "unknown" ;constrain.bbp= F
# constrain.shallow.bgc = F ; constrain.shallow.iop = T; manual_par0 = T
# #type_Rrs_below = "shallow"
# 
# inv_bound = create_init_bound(rrs_inv_type = type_Rrs_below, manual_par0 = manual_par0, 
#                               constrain.bbp = constrain.bbp, 
#                               constrain.shallow.bgc = constrain.shallow.bgc, 
#                               constrain.shallow.iop = constrain.shallow.iop, 
#                               pop.sd =  pop.sd, 
#                               init_par = c(sapply(param_vec, mean ), 0.5,0.5,0.5,0.05),
#                               upper_par = c(sapply(param_vec, max), 1,1,1,10),
#                               lower_par = c(sapply(param_vec, min ), 0.0,0.0,0.0,0.00001)
# )
# 
# par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
# lower.bound = inv_bound$lower.bound
# 
# doOptimization_deep_const(obsdata = obsdata, 
#                           bbp_const_val = 0.007,
#                             max_par = upper.bound,
#                             min_par = lower.bound, 
#                             init_par = par0, wavelength_sim = wavelength, 
#                             QAA_mode = F,
#                             qaa_prefit = F, qaa_slope = F, manual_slope = T, 
#                             manual_slope_vals = c("s_g"= 0.014, "s_d"= 0.003, "gamma"=0.5),
#                             opt_method = methods.opt[4],
#                             sa_model = "am03", 
#                             obj_fn = obj.fn )
# 
# doOptimization_deep_unconst(obsdata = obsdata, 
#                             max_par = upper.bound,
#                             min_par = lower.bound, 
#                             init_par = par0, wavelength_sim = wavelength, 
#                             QAA_mode = F,
#                             qaa_prefit = F, qaa_slope = F, manual_slope = T, 
#                             manual_slope_vals = c("s_g"= 0.014, "s_d"= 0.003, "gamma"=0.5),
#                             opt_method = methods.opt[4],
#                             sa_model = "am03", 
#                             obj_fn = obj.fn )
# 
# doOptimization_shallow_BGC_const(obsdata = obsdata, 
#                             max_par = upper.bound,
#                             min_par = lower.bound, bgc_const_val = c(4.96, 1, 0.007),
#                             init_par = par0, wavelength_sim = wavelength, 
#                             QAA_mode = F,
#                             qaa_prefit = F, qaa_slope = F, manual_slope = F, 
#                             manual_slope_vals = c(0.014, 0.003, 0.5),
#                             opt_method = methods.opt[4],
#                             sa_model = "am03", 
#                             obj_fn = obj.fn )
# 
# doOptimization_shallow_IOP_const(obsdata = obsdata, 
#                                  max_par = upper.bound,
#                                  min_par = lower.bound, bgc_const_val = c(4.96, 1, 0.007),
#                                  iop_const_path = config_list$` constrain_iop`,
#                                  init_par = par0, wavelength_sim = wavelength, 
#                                  QAA_mode = F,
#                                  qaa_prefit = F, qaa_slope = F, manual_slope = F, 
#                                  manual_slope_vals = c(0.014, 0.003, 0.5),
#                                  opt_method = methods.opt[4],
#                                  sa_model = "am03", 
#                                  obj_fn = obj.fn )
# 
# doOptimization_shallow_unconst(obsdata = obsdata, 
#                             max_par = upper.bound,
#                             min_par = lower.bound, bgc_const_val = c(4.96, 1, 0.007),
#                             init_par = par0, wavelength_sim = wavelength, 
#                             QAA_mode = F,
#                             qaa_prefit = F, qaa_slope = F, manual_slope = F, 
#                             manual_slope_vals = c(0.014, 0.003, 0.5),
#                             opt_method = methods.opt[4],
#                             sa_model = "am03", 
#                             obj_fn = obj.fn )

########################################################################################
inverse_runGrad <- function(obsdata, rrs_type,
                            max_par, min_par, init_par,
                            param_lab = c("chl","adg443","bbp555", 
                                          "H",  rep(paste0("fa", (1:rb_count))), 
                                          "pop_sd"),
                            
                            constrain_config = c("const_bbp" = F, "cost_bgc"=F, "cost_IOP"=F),
                            bgc_const_val, iop_const_path, bbp_const_val,
                            
                            qaa_prefit, QAA_mode,
                            
                            qaa_slope, manual_slope, 
                            manual_slope_vals= c("s_g"=0.014, "s_d"=0.003, "gamma"=0.5),
                            wavelength_sim, 
                            sa_model, obj_fn, opt_method,
                            plot_rrs){
  
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### GET THAT GRADIENT #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  cat(paste0("\033[0;34m MODEL CONSTRAINTS ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Backscatter constrain: ",constrain_config[1],"\033[0m","\n"))
  if (constrain_config[1] == T) {
    cat(paste0("\033[0;32m","bbp(555) constrained value: ",bbp_const_val,"\033[0m","\n"))
  }
  
  cat(paste0("\033[0;32m","Biogeochemical variable constrain: ",constrain_config[2],"\033[0m","\n"))
  if (constrain_config[2] == T) {
    cat(paste0("\033[0;32m","BGC constrained value: ",bgc_const_val,"\033[0m","\n"))
  }
  cat(paste0("\033[0;32m","IOP constrain: ",constrain_config[3],"\033[0m","\n"))
  if (constrain_config[3] == T) {
    cat(paste0("\033[0;32m","Spectral IOP paths: ",paste0(iop_const_path, collapse = ","),"\033[0m","\n"))
  }
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m SPECTRAL SLOPE SETTINGS ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope for QAA :=",qaa_slope,"\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope from user:=",manual_slope,"\033[0m","\n"))
  
  if (manual_slope == T) {
    cat(paste0("\033[0;32m","Slope values:=",manual_slope_vals,"\033[0m","\n"))
    
  }
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m ANCILIARY ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Rrs type: ",rrs_type,"\033[0m","\n"))
  
  if (sa_model == "am03") {
    cat(paste0("\033[0;32m","Paramterization for forward model: Albert & Mobley (2003) \033[0m","\n"))
    
  } else {
    cat(paste0("\033[0;32m","Paramterization for forward model: Lee et al., (1998) \033[0m","\n"))
    
  }
  cat(paste0("\033[0;32m","QAA Mode: ",QAA_mode,"\033[0m","\n"))
  cat(paste0("\033[0;32m","QAA Prefit:",qaa_prefit,"\033[0m","\n"))
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  Sys.sleep(1)
  
  
  
  inv_bound = create_init_bound(rrs_inv_type = rrs_type, manual_par0 = T, 
                                constrain.bbp = constrain_config[1], 
                                constrain.shallow.bgc = constrain_config[2], 
                                constrain.shallow.iop = constrain_config[3], 
                                pop.sd =  "unknown", 
                                init_par = init_par,
                                upper_par = max_par,
                                lower_par = min_par
  )
  
  init_par = inv_bound$par0 ; max_par = inv_bound$upper.bound  
  min_par = inv_bound$lower.bound
  
  #browser()
  
  if (rrs_type == "deep") {
    type_Rrs_below = "deep"
    
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    #browser()
    
    if (constrain_config[1] == TRUE && (all(constrain_config[-1])) == FALSE) {
      
      inverse_output = doOptimization_deep_const(obsdata = obsdata,
                                                   max_par = max_par,
                                                   min_par = min_par, 
                                                   init_par = init_par,
                                                   wavelength_sim = wavelength_sim, 
                                                   QAA_mode = QAA_mode, 
                                                  bbp_const_val = bbp_const_val,
                                                   qaa_prefit = qaa_prefit, 
                                                   qaa_slope = qaa_slope, 
                                                   manual_slope = manual_slope, 
                                                   manual_slope_vals = manual_slope_vals,
                                                   opt_method = opt_method,
                                                   sa_model = sa_model, 
                                                   obj_fn = obj.fn )
      
    } else {
      
      inverse_output = doOptimization_deep_unconst(obsdata = obsdata,
                                                   max_par = max_par,
                                                   min_par = min_par, 
                                                   init_par = init_par, wavelength_sim = wavelength_sim, 
                                                   QAA_mode = QAA_mode,
                                                   qaa_prefit = qaa_prefit, 
                                                   qaa_slope = qaa_slope, manual_slope = manual_slope, 
                                                   manual_slope_vals = manual_slope_vals,
                                                   opt_method = opt_method,
                                                   sa_model = sa_model, 
                                                   obj_fn = obj.fn )
      
    }
    
    }
    
  
  if (rrs_type == "shallow") {
    
    type_Rrs_below = "shallow"
    
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    
    if (constrain_config[3] == TRUE && (all(constrain_config[-3])) == FALSE) {
      
      inverse_output = doOptimization_shallow_IOP_const(obsdata = obsdata, 
                                            max_par = max_par, 
                                            iop_const_path = iop_const_path,
                                            min_par = min_par, bgc_const_val = bgc_const_val,
                                            init_par = init_par, wavelength_sim = wavelength_sim, 
                                            QAA_mode = QAA_mode,
                                            qaa_prefit = qaa_prefit, 
                                            qaa_slope = qaa_slope, manual_slope = manual_slope, 
                                            manual_slope_vals = manual_slope_vals,
                                            opt_method = opt_method,
                                            sa_model = sa_model, 
                                            obj_fn = obj.fn)
      
    } else {
      
      if (constrain_config[2] == TRUE && (all(constrain_config[-2])) == FALSE) {
        
        inverse_output = doOptimization_shallow_BGC_const(obsdata = obsdata, 
                                                          max_par = max_par,
                                                          min_par = min_par, bgc_const_val = bgc_const_val,
                                                          init_par = init_par, wavelength_sim = wavelength_sim, 
                                                          QAA_mode = QAA_mode,
                                                          qaa_prefit = qaa_prefit, 
                                                          qaa_slope = qaa_slope, manual_slope = manual_slope, 
                                                          manual_slope_vals = manual_slope_vals,
                                                          opt_method = opt_method,
                                                          sa_model = sa_model, 
                                                          obj_fn = obj.fn)
        
      } else {
        
        if (all(constrain_config) == F) {
          
          inverse_output = doOptimization_shallow_unconst(obsdata = obsdata, 
                                                          max_par = max_par,
                                                          min_par = min_par, bgc_const_val = bgc_const_val,
                                                          init_par = init_par, wavelength_sim = wavelength_sim, 
                                                          QAA_mode = QAA_mode,
                                                          qaa_prefit = qaa_prefit, 
                                                          qaa_slope = qaa_slope, manual_slope = manual_slope, 
                                                          manual_slope_vals = manual_slope_vals,
                                                          opt_method = opt_method,
                                                          sa_model = sa_model, 
                                                          obj_fn = obj.fn  )
          
        }
        
      }
      
    }
    
  }
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### INVERSION IS COMPLETED #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  names(inverse_output) = c(param_lab, paste0("sd_", param_lab))
  return(inverse_output)
}

# #Instantiate inversion objective function and the optimization scheme 
# pop.sd = "unknown" ;constrain.bbp= F
# constrain.shallow.bgc = F ; constrain.shallow.iop = F; manual_par0 = T
# type_Rrs_below = "shallow"
# 
# inv_bound = create_init_bound(rrs_inv_type = type_Rrs_below, manual_par0 = manual_par0, 
#                               constrain.bbp = constrain.bbp, 
#                               constrain.shallow.bgc = constrain.shallow.bgc, 
#                               constrain.shallow.iop = constrain.shallow.iop, 
#                               pop.sd =  pop.sd, 
#                               init_par = c(2,1,0.007,5,0.5,0.5,0.5,0.05),
#                               upper_par = c(sapply(param_vec, max), 1,1,1,10),
#                               lower_par = c(sapply(param_vec, min ),0.0,0.0,0.0,0.00001)
# )
# 
# par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
# lower.bound = inv_bound$lower.bound
# 
# test_inv = inverse_runGrad(obsdata = obsdata, rrs_type = type_Rrs_below,
#                 max_par = upper.bound,
#                 min_par = lower.bound,
#                 init_par = par0,
#                 param_lab = c("chl","adg443","bbp555",
#                               "H",  rep(paste0("fa", (1:rb_count))),
#                               "pop_sd"),
#                 constrain_config = c(F, F, F), bgc_const_val = c(4.96,1.003,0.007),
#                 iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
#                                    "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
#                 qaa_prefit = F, QAA_mode = F, qaa_slope = F, manual_slope = F,
#                 manual_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=0.46),
#                 wavelength_sim =  wavelength,
#                 sa_model = "am03", obj_fn = obj.fn, opt_method = methods.opt[4]
#                 )

###############################################################################################################

#===================================================================
# Unconstrained inversion for deep water (with QAA support)
#===================================================================
doOptimization_deep_unconst_back <- function(obsdata, par0_base, wl, qaa_prefit, QAA_mode, 
                                        sa_model, obj_fn, sub_surf_translate) {
  
  #Assign the inversion vector of rrs
  rrs_inverse_input = obsdata
  
  if (QAA_mode == F) {
    
    if (qaa_prefit == TRUE) {
      print("QAA PREFIT : ON")
      qaa_output = QAA.v5(waves = wavelength, Rrs = rrs_inverse_input, 
                          translate_subsurf = sub_surf_translate)
      par0 = c(chl = qaa_output$chl, adg440 =qaa_output$a_dg_443, 
               bbp550 = qaa_output$b_bp_555, "pop_sd" = par0_base[length(par0_base)])
      
      if (any(par0 == Inf) | any(is.na(par0))) {
        print("Prefit failed, default values wil be used")
        par0 = par0_base
        
        upper.bound <- upper.bound_base  
        lower.bound <- lower.bound_base
      } else {
        print(paste0("prefit values are computed as chl: ", signif(par0[1], digits = 3), 
                     " adg443: ", signif(par0[2], digits = 3), 
                     " bbp555: ",signif(par0[3], digits = 3), 
                     " pop_sd: ",signif(par0[4], digits = 3)))
        
        upper.bound <- c(par0[1:2] + 0.15*par0[1:2], upper.bound_base[3:4])  
        lower.bound <- c(par0[1:2] - 0.15*par0[1:2], lower.bound_base[3:4])
      }
      
    } else {
      print("QAA PREFIT : OFF")
      par0 = par0_base
      print(paste0("prefit values are computed as chl: ", signif(par0[1], digits = 3), 
                   " adg443: ", signif(par0[2], digits = 3), 
                   " bbp555: ",signif(par0[3], digits = 3), 
                   " pop_sd: ",signif(par0[4], digits = 3)))
      upper.bound <- upper.bound_base  
      lower.bound <- lower.bound_base
    }
    
    inverse_output_deep <- suppressWarnings(solve.objective.inverse.deep.final.fast(
      wave = wl,
      initial = par0, 
      bbp.constrain  = F,
      bbp.constrain.value =  const.input[1],
      
      obsdata = as.numeric(rrs_inverse_input),
      
      sa.model = sa_model, 
      obj.fn = obj_fn,
      
      auto_spectral_slope = T,
      manual_spectral_slope = F, 
      
      manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, 
                                     "gamma"=refexponent),
      
      method.opt = methods.opt[4],
      lower.b = lower.bound,
      upper.b = upper.bound 
      #,batch = FALSE, pop.sd = FALSE
    ))
    
    Fit.optimized.ssobj.batch <- c(inverse_output_deep[[1]]$estimates,
                                   inverse_output_deep[[1]]$`sd(+/-)`)
    
  } else {
    
    
    qaa_output = QAA.v5(waves = wl, Rrs = rrs_inverse_input, 
                        translate_subsurf = sub_surf_translate)
    Fit.optimized.ssobj.batch = c(chl = qaa_output$chl, adg443 =qaa_output$a_dg_443, 
             bbp555 = qaa_output$b_bp_555)
    
    if (any(Fit.optimized.ssobj.batch == Inf) | any(is.na(Fit.optimized.ssobj.batch))) {
      print("QAA failed, default values wil be used")
      Fit.optimized.ssobj.batch = par0_base[-length(par0_base)]
  }
  
  
  return(Fit.optimized.ssobj.batch)
  }
}

# doOptimization_deep_unconst(obsdata = obsdata, par0_base = par0, wl = wavelength, 
#                             qaa_prefit = F, QAA_mode = F, sa_model = "am03", obj_fn = obj.fn)

#===================================================================
### Unconstrained inversion for shallow water
#===================================================================
doOptimization_shallow_unconst_back <- function(obsdata, par0, wl, sa_model, obj_fn,opt_method) {
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
    
    wave = wl,
    
    constrain.shallow.iop = F,  
    
    unconstrained = T,
    
    constrained.bgc = F, 
    constrain.bgc.value = c(Fit.input$chl, Fit.input$acdom.440+ 
                              Fit.input$anap.440),
    
    
    initial = as.numeric(par0), 
    obsdata = as.numeric(obsdata),
    
    auto_spectral_slope = F,
    manual_spectral_slope = F, 
    
    manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.01160, "gamma"=0.46),
    
    sa.model = sa_model, 

    obj.fn = obj_fn,
    method.opt = opt_method,
    
    lower.b = lower.bound,
    upper.b = upper.bound, 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}

#===================================================================
# BGC constrained inversion for shallow water
#===================================================================
doOptimization_shallow_BGC_const_back <- function(obsdata, par0, wl,sa_model, obj_fn, opt_method) {
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
    
    wave = wl,
    
    unconstrained = F, 
    constrain.shallow.iop = F, abs_path = NA, bb_path = NA,
    constrained.bgc  = T, 
    
    constrain.bgc.value = as.numeric(obsdata[(length(obsdata)-3): (length(obsdata)-1)]),
    
    initial = as.numeric(par0), 
    obsdata = as.numeric(obsdata[1:(length(obsdata)-4)]),
    
    auto_spectral_slope = F,
    manual_spectral_slope = F, 
    
    manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=1),
    
    sa.model = sa_model, 

    obj.fn = obj_fn,
    method.opt = opt_method,
    
    lower.b = lower.bound,
    upper.b = upper.bound, 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}

# inverse_parallel_shallow <- function(constr = F, input_list, par0) {
#   
#   load("./global_var.RData")
#   
#   # Unconstrained inversion for shallow water
#   doOptimization_shallow_unconst <- function(obsdata, par0) {
#     inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
#       constrain.shallow.iop = F, 
#       
#       unconstrained = T,
#       
#       constrained.bgc = F, 
#       constrain.bgc.value = c(Fit.input$chl, Fit.input$acdom.440+ 
#                                 Fit.input$anap.440),
#       
#       
#       initial = as.numeric(par0), 
#       obsdata = as.numeric(obsdata),
#       
#       auto_spectral_slope = F,
#       manual_spectral_slope = F, 
#       
#       manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=1),
#       
#       sa.model = "am03", 
#       
#       obj.fn = obj[1],
#       method.opt = methods.opt[4],
#       
#       lower.b = lower.bound,
#       upper.b = upper.bound, 
#       wave = wavelength
#     ))
#     
#     Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
#                                    inverse_output[[1]]$`sd(+/-)`)
#     return(Fit.optimized.ssobj.batch)
#   }
#   
#   # BGC constrained inversion for shallow water
#   doOptimization_shallow_BGC_const <- function(obsdata, par0) {
#     inverse_output <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
#       unconstrained = F, 
#       constrain.shallow.iop = F, abs_path = NA, bb_path = NA,
#       constrained.bgc  = T, 
#       constrain.bgc.value = as.numeric(obsdata[(length(obsdata)-3): (length(obsdata)-1)]),
#       
#       initial = as.numeric(par0), 
#       obsdata = as.numeric(obsdata[1:(length(obsdata)-4)]),
#       
#       auto_spectral_slope = F,
#       manual_spectral_slope = F, 
#       
#       manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=1),
#       
#       sa.model = "am03", 
#       
#       obj.fn = obj[1],
#       method.opt = methods.opt[4],
#       
#       lower.b = lower.bound,
#       upper.b = upper.bound, 
#       wave = wavelength
#     ))
#     
#     Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
#                                    inverse_output[[1]]$`sd(+/-)`)
#     return(Fit.optimized.ssobj.batch)
#   }
#   
#   # Set up parallel processing
#   numCores <- detectCores()
#   print(paste0("Number of system processor found: ", numCores))
#   
#   cl <- makeCluster(detectCores() - 1)
#   registerDoParallel(cl)
#   print(paste0("The Cluster is created with ", length(cl), " number of processors"))
#   
#   start_time = Sys.time()
#   if (constr == TRUE) {
#     
#     cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
#     cat(paste0("\033[0;39m","########### CONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
#     cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
#     # Apply the function to each element in the lists using %dopar%
#     output_list <- foreach(i = names(input_list)
#                            ,
#                            .export = c("obj", "methods.opt","lower.bound", 
#                                        "upper.bound", "wavelength",
#                                        "type_case_water", "type_Rrs_below", "type_Rrs_water"),
#                            .verbose = T,
#                            #.combine = rbind,
#                            .packages="data.table"
#     ) %dopar% {
#       source("./R/SABER_forward_fast.R")
#       source("./R/solve.objective.inverse_fast.R")
#       source("./R/saber.inversion.vector.parallel.R")
#       apply(input_list[[i]], 1, doOptimization_shallow_BGC_const,
#             par0 = par0
#             
#       )
#     }
#     
#   } else {
#     cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
#     cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
#     cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
#     # Apply the function to each element in the lists using %dopar%
#     output_list <- foreach(i = names(input_list)
#                            ,  
#                            .export = c("obj", "methods.opt","lower.bound", 
#                                        "upper.bound", "wavelength",
#                                        "type_case_water", "type_Rrs_below", "type_Rrs_water"),
#                            #.combine = rbind,
#                            .verbose = T,
#                            .packages="data.table"
#     ) %dopar% {
#       source("./R/SABER_forward_fast.R")
#       source("./R/solve.objective.inverse_fast.R")
#       source("./R/saber.inversion.vector.parallel.R")
#       apply(input_list[[i]], 1, doOptimization_shallow_unconst,
#             par0 = par0
#             
#       )
#     }
#     
#   }
#   end.time = Sys.time(); time_taken <- end.time - start.time
#   cat(paste0("\033[0;33m","Time Taken for the inversion: ",time_taken," secs. \033[0m","\n"))
#   return(output_list)
#   stopCluster(cl)
# }
