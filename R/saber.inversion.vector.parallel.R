#-----------------------------------------------------------------------------------------
#Vectorize the inversion for multiple Rrs observation
#-----------------------------------------------------------------------------------------
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

# Unconstrained inversion for shallow water
doOptimization_shallow_unconst <- function(obsdata, par0, wl) {
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
    manual_spectral_slope = T, 
    
    manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.01160, "gamma"=0.46),
    
    sa.model = "am03", 

    obj.fn = obj[1],
    method.opt = methods.opt[4],
    
    lower.b = lower.bound,
    upper.b = upper.bound, 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}

# BGC constrained inversion for shallow water
doOptimization_shallow_BGC_const <- function(obsdata, par0, wl) {
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
    
    sa.model = "am03", 

    obj.fn = obj[1],
    method.opt = methods.opt[4],
    
    lower.b = lower.bound,
    upper.b = upper.bound, 
  ))
  
  Fit.optimized.ssobj.batch <- c(inverse_output[[1]]$estimates,
                                 inverse_output[[1]]$`sd(+/-)`)
  return(Fit.optimized.ssobj.batch)
}
