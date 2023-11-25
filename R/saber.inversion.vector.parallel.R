#-----------------------------------------------------------------------------------------
#Vectorize the inversion for multiple Rrs observation
#-----------------------------------------------------------------------------------------

#===================================================================
### Deep water unconstrained inversion
#===================================================================


# Unconstrained inversion for deep water (with QAA support)
doOptimization_deep_unconst <- function(obsdata, par0_base, wl, qaa_prefit, QAA_mode, 
                                        sa_model, obj_fn) {
  
  #Assign the inversion vector of rrs
  rrs_inverse_input = obsdata
  
  if (QAA_mode == F) {
    
    if (qaa_prefit == TRUE) {
      print("QAA PREFIT : ON")
      qaa_output = QAA.v5(waves = wavelength, Rrs = rrs_inverse_input)
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
    
    
    qaa_output = QAA.v5(waves = wl, Rrs = rrs_inverse_input)
    Fit.optimized.ssobj.batch = c(chl = qaa_output$chl, adg443 =qaa_output$a_dg_443, 
             bbp555 = qaa_output$b_bp_555)
    
    if (any(Fit.optimized.ssobj.batch == Inf) | any(is.na(Fit.optimized.ssobj.batch))) {
      print("QAA failed, default values wil be used")
      Fit.optimized.ssobj.batch = par0_base[-length(par0_base)]
  }
  
  
  return(Fit.optimized.ssobj.batch)
  }
}

#===================================================================
### Unconstrained inversion for shallow water
#===================================================================
doOptimization_shallow_unconst <- function(obsdata, par0, wl, sa_model, obj_fn,opt_method) {
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
doOptimization_shallow_BGC_const <- function(obsdata, par0, wl,sa_model, obj_fn, opt_method) {
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
