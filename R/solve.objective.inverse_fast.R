#===========================================================================
# solve.objective.inverse.final.R solves the inverse objective function  to retrieve the water components 
#along with bathymetry and bottom reflectance (for shallow water) given the initial param vals and Rrs.

#There are three optimization methods:

#1. Minimize the residual sum-square-error using non-linear multivariate optimization using 
#non-gradient based method
#2. Maximize the conditional log-likelihood (P(Data|Theta) of the model noise and obtain hessian
#derivative of the function, iff the function is differentiable over paramteric space)
#3. Spectral Error Index from Lee et al. 1999 

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================

#===============================================================================================================
#Solver for shallow waters
#===============================================================================================================

#Formula to implement SAM objective function
# acos((sum(sum(rrs.forward.am.param.conc.dg_comp_sicf_fdom*rrs.forward.am)))/
#        ((sum(rrs.forward.am.param.conc.dg_comp_sicf_fdom^2)^0.5)*(sum(rrs.forward.am^2)^0.5)))

solve.objective.inverse.shallow.final.fast <- function(
    constrain.shallow.iop,
    abs_path,
    bb_path,
    
    
    constrained.bgc,
    constrain.bgc.param = c("chl", "adg", "bbp"), 
    constrain.bgc.value,
    
    unconstrained,
    
    auto_spectral_slope,
    manual_spectral_slope,
    manual_spectral_slope_vals = c("s_g"=0.015, "s_d"=0.01160, "gamma"=0.5),
    
    obj.fn, initial, sa.model="am03", method.opt,
    
    obsdata,  wave = wavelength,
    
    lower.b, upper.b
     
    #,batch=FALSE, pop.sd=FALSE
    ){
  
  #initial_rb_length = length(initial[5:(length(initial)-1)])
  initial_rb_length = 3
  
  
  if (auto_spectral_slope == TRUE & manual_spectral_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  myFun <- function(x) {
    NA
  }
  
  #==========================================================================
  # For IOP constrained inversion
  #==========================================================================
  if (constrain.shallow.iop == TRUE) {
    initial_rb_length = length(initial[2:(length(initial)-1)])
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### IOP CONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    Sys.sleep(1)
    if (obj.fn == "log-LL") {
      print("Log-Likelihood will be used to construct objective function")
      ##Create log-likelihood function
      
      NLL_unconstr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = T, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            #chl = pars[1], 
            #a_dg = pars[2],
            #bbp.550 = pars[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave,plot = F
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = T, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            #chl = pars[1], 
            #a_dg = pars[2],
            #bbp.550 = pars[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave, plot = F
          )
        }
        
        # Negative log-likelihood
        smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                           log = TRUE))
        return(smull)
      }
    } else {
      
      print("Spectral Error index from Lee et al. 1999 will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_unconstr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
            Gpred = Saber_forward_fast(
              use_true_IOPs = T, 
              a_non_water_path = abs_path,
              bb_non_water_path = bb_path,
              
              #chl = pars[1], 
              #a_dg = pars[2],
              #bbp.550 = pars[3],
              
              z = pars[1],
              rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
              
              
              slope.parametric = auto_spectral_slope,
              Rrs_input_for_slope = data, 
              
              use_manual_slope =manual_spectral_slope,
              manual_slope =  manual_spectral_slope_vals,
              
              verbose = F, wavelength = wave, plot = F
            )
            
          
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = T, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            #chl = pars[1], 
            #a_dg = pars[2],
            #bbp.550 = pars[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave, plot = F
          )
        }
        
        rrs_est = Gpred[[1]]$Rrs
        
        # Sum-squared residual of error (SSE)
        #sse= sum((data - rrs_est)^2)
        #return(sse)
        
        #The Spectral error index from Lee 1999
        # Define the spectral regions
        region1 <- which(wavelength >= 400 & wavelength <= 675)
        region2 <- which(wavelength >= 750 & wavelength <= 830)
        
        # Calculate the numerator of the error index
        numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                            sum((data[region2] - rrs_est[region2])^2))
        
        # Calculate the denominator of the error index
        denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
        
        # Calculate the error index
        smull <- numerator / denominator
        
        return(smull)
        
      }
      
    }
    
    ##Optimize the error function
    cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
    
    if (auto_spectral_slope == TRUE) {
      print("Spectral slopes will be calculated using QAA")
    } else {
      if (manual_spectral_slope == TRUE) {
        print(paste0("Spectral slopes will be used from user supplied values: ",
                     manual_spectral_slope_vals[1], ", " , manual_spectral_slope_vals[2], ", " ,manual_spectral_slope_vals[3]))
      } else {
        print("Spectral slopes will be used same as model default")
      }
    }
    
    #Instantiate initial values
    par0 = c(z = initial[1], 
             initial[2:(1+initial_rb_length)],
             pop.sd = initial[length(initial)])
    
    
    
    #cat(paste0("\033[0;33m","####################CONSTRAINED INVERSION#########################","\033[0m","\n"))
    
    cat(paste0("\033[0;32m","Initial values are: zB=", par0[1],",  RB={", toString(as.numeric(par0[2:(length(par0)-1)])),"},  population.sigma=", par0[length(par0)],"\033[0m","\n"))
    
    
    cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
    #Sys.sleep(2)
    start.time = Sys.time()
    
    if (sa.model == "am03") {
      print("Albert & Mobley 2003 SA model used for forward model")
    } else {
      print("Lee et al. 1999 SA model used for forward model")
    }
    
    if (method.opt == "L-BFGS-B") {
      
      print("L-BFGS-B Optimization will be used for inversion")
      MLE_estimates = optim(par = as.numeric(par0), fn = NLL_unconstr, data = obsdata, 
                            lower = as.numeric(lower.b),     # Lower bound on parameters
                            upper = as.numeric(upper.b),  # Upper bound on parameters
                            #gr = grad.calculate,
                            method = method.opt,
                            control = list(#fnscale = 1, 
                              parscale = as.numeric(abs(par0))),
                            #hessian = TRUE
      )
    } 
    
    if (method.opt == "levenberg-marqardt") {
      
      print("Levenberg-Marquardt Optimization will be used for inversion")
      LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_unconstr,data = obsdata, print.info = F)
      MLE_estimates = data.frame("par"=LM_estimates$b)
    }
    
    if (method.opt == "auglag") {
      print("Augmented Lagriangian with equality constraints will be used for inversion")
      
      
      fheq <- function(pars,data) sum(pars[2:(length(pars)-1)]) - 1
      
      fhin <- function(pars,data) c(
        #bounds for aerial fractions
        pars[2:(length(pars)-1)]
      )
      
      MLE_estimates <- alabama::auglag(fn = NLL_unconstr,par = par0, 
                                       #gr = pracma::grad(NLL_unconstr, x0 = par, data= obsdata),
                                       heq = fheq,
                                       hin = fhin,
                                       #lower = lower.bound,
                                       #upper = upper.bound,
                                       #data=surface_rrs_translate(insitu.data),
                                       data = obsdata,
                                       control.outer = list(trace = F, method = "nlminb")
      )
      
      
      
      
      
      #print(inverse_output_fmincon$par, digits=3)
      
      #print(sum(inverse_output_fmincon$par[5:9]))
      
      #sprintf("%.5f", inverse_output_fmincon$par)
      
    }  
    
    if (method.opt == "Nelder-Mead" |method.opt ==  "SANN" | method.opt ==  "Brent") {
      
      print("Simplex based Optimization will be used for inversion")
      MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata, 
                            #lower = c(0, 0, 0),     # Lower bound on parameters
                            #upper = c(10, 2, 0.3),  # Upper bound on parameters
                            method = method.opt,
                            control = list(#fnscale = 1, 
                              parscale = as.numeric(abs(par0_constrained))),
                            hessian = FALSE)
      
      
    }
    
    
    cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
    
    #print(MLE_estimates$par)
    
    cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
    Sys.sleep(1)
    
    #Calculate hessian matrix for var-covar matrix
    if (method.opt == "auglag") {
      hessian.inverse <- MLE_estimates$hessian 
    } else {
      hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_unconstr, data=obsdata)
    }
    
    rownames(hessian.inverse) <- names(par0)
    colnames(hessian.inverse) <- names(par0)
    cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
    prmatrix(hessian.inverse)
    cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
    
    MLE_par <- MLE_estimates$par
    
    result <- tryCatch({
      sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
      error = myFun
    )
    if (!is.numeric(result)) {
      cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
      MLE_SE = result
    } else {
      MLE_SE <- result
    }
    MLE <- data.table("param" = names(par0),
                      "estimates" = MLE_par,
                      "sd(+/-)" = MLE_SE)
    
    print("The retrieved parameters are:")
    prmatrix(MLE)
    
    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
      cat(paste0("\033[0;32m","#################### CONVERGENCE: GLOBAL #####################","\033[0m","\n"))
    } else {
      convergence = "FALSE"
      cat(paste0("\033[0;34m","#################### CONVERGENCE: LOCAL #####################","\033[0m","\n"))
    }
    end.time = Sys.time(); time_taken <- end.time - start.time
    return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
    
    cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
  }
  
  #==========================================================================
  # For BGC variable constrained inversion
  #==========================================================================
  if (constrained.bgc == TRUE) {
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### BGC CONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    Sys.sleep(1)
    
    #initial <- abs(initial)
    if (obj.fn == "log-LL") {
      print("Log-Likelihood will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_constr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = constrain.bgc.value[1], 
            a_dg = constrain.bgc.value[2],
            bbp.550 = constrain.bgc.value[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = constrain.bgc.value[1], 
            a_dg = constrain.bgc.value[2],
            bbp.550 = constrain.bgc.value[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        # Negative log-likelihood
        
        smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                               log = TRUE)) 
        return(smull)
      }
    } else {
      print("Spectral Error index from Lee et al. 1999 will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_constr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = constrain.bgc.value[1], 
            a_dg = constrain.bgc.value[2],
            bbp.550 = constrain.bgc.value[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = constrain.bgc.value[1], 
            a_dg = constrain.bgc.value[2],
            bbp.550 = constrain.bgc.value[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        rrs_est = Gpred[[1]]$Rrs
        
        # Sum-squared residual of error (SSE)
        #sse= sum((data - rrs_est)^2)
        #return(sse)
        
        #The Spectral error index from Lee 1999
        # Define the spectral regions
        region1 <- which(wavelength >= 400 & wavelength <= 675)
        region2 <- which(wavelength >= 750 & wavelength <= 830)
        
        # Calculate the numerator of the error index
        numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                            sum((data[region2] - rrs_est[region2])^2))
        
        # Calculate the denominator of the error index
        denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
        
        # Calculate the error index
        smull <- numerator / denominator
        
        return(smull)
        
      }
      
    }
      
      ##Optimize the error function
      cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
      
      if (auto_spectral_slope == TRUE) {
        print("Spectral slopes will be calculated using QAA")
      } else {
        if (manual_spectral_slope == TRUE) {
          print(paste0("Spectral slopes will be used from user supplied values: ",
                       manual_spectral_slope_vals[1], ", " , manual_spectral_slope_vals[2], ", " ,manual_spectral_slope_vals[3]))
        } else {
          print("Spectral slopes will be used same as model default")
        }
      }
      
      #Instantiate initial values
      par0 = c(#chl = initial[1], adg440 = initial[2], bbp550 = initial[3], 
                 
                z = initial[1], 
                initial[2:(1+initial_rb_length)],
                pop.sd = initial[length(initial)])
      
      par0_constrained = par0
      
      
      cat(paste0("\033[0;33m","####################CONSTRAINED INVERSION#########################","\033[0m","\n"))
      
      cat(paste0("\033[0;32m","Initial values are: zB=", par0[1],",  RB={", toString(as.numeric(par0[2:(length(par0)-1)])),"},  population.sigma=", par0[length(par0)],"\033[0m","\n"))
      
      
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      #Sys.sleep(2)
      start.time = Sys.time()
      
      if (sa.model == "am03") {
        print("Albert & Mobley 2003 SA model used for forward model")
      } else {
        print("Lee et al. 1999 SA model used for forward model")
      }
      
      if (method.opt == "L-BFGS-B") {
        
        print("L-BFGS-B Optimization will be used for inversion")
        MLE_estimates = optim(par = as.numeric(par0_constrained), fn = NLL_constr, data = obsdata, 
                              lower = as.numeric(lower.b),     # Lower bound on parameters
                              upper = as.numeric(upper.b),  # Upper bound on parameters
                              #gr = grad.calculate,
                              method = method.opt,
                              control = list(#fnscale = 1, 
                                parscale = as.numeric(abs(par0_constrained))),
                              #hessian = TRUE
        )
      } 
      
      if (method.opt == "levenberg-marqardt") {
        
        print("Levenberg-Marquardt Optimization will be used for inversion")
        LM = marqLevAlg::marqLevAlg(b = par0_constrained, fn = NLL_constr,data = obsdata, print.info = F)
        MLE_estimates = data.frame("par"=LM_estimates$b)
      }
       
      if (method.opt == "auglag") {
        print("Augmented Lagriangian with equality constraints will be used for inversion")
        
        
        fheq <- function(pars,data) sum(pars[2:(length(pars)-1)]) - 1
        
        fhin <- function(pars,data) c(
          #bounds for aerial fractions
          pars[2:(length(pars)-1)]
        )
        
        MLE_estimates <- alabama::auglag(fn = NLL_constr,par = par0_constrained, 
                                         #gr = pracma::grad(NLL_unconstr, x0 = par, data= obsdata),
                                         heq = fheq,
                                         hin = fhin,
                                         #lower = lower.bound,
                                         #upper = upper.bound,
                                         #data=surface_rrs_translate(insitu.data),
                                         data = obsdata,
                                         control.outer = list(trace = F, method = "nlminb")
        )
        
        
        
        
        
        #print(inverse_output_fmincon$par, digits=3)
        
        #print(sum(inverse_output_fmincon$par[5:9]))
        
        #sprintf("%.5f", inverse_output_fmincon$par)
        
      }  
      
      if (method.opt == "Nelder-Mead" |method.opt ==  "SANN" | method.opt ==  "Brent") {
        
        print("Simplex based Optimization will be used for inversion")
        MLE_estimates = optim(par = par0_constrained, fn = NLL_constr, data = obsdata, 
                              #lower = c(0, 0, 0),     # Lower bound on parameters
                              #upper = c(10, 2, 0.3),  # Upper bound on parameters
                              method = method.opt,
                              control = list(#fnscale = 1, 
                                parscale = as.numeric(abs(par0_constrained))),
                              hessian = FALSE)
        
        
      }
      
      
      cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
      
      #print(MLE_estimates$par)
      
      cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
      Sys.sleep(1)
      
      #Calculate hessian matrix for var-covar matrix
      if (method.opt == "auglag") {
        hessian.inverse <- MLE_estimates$hessian 
      } else {
        hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_constr, data=obsdata)
      }
      
      rownames(hessian.inverse) <- names(par0_constrained)
      colnames(hessian.inverse) <- names(par0_constrained)
      cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
      prmatrix(hessian.inverse)
      cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
      
      MLE_par <- MLE_estimates$par
      
      result <- tryCatch({
        sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
        error = myFun
      )
      if (!is.numeric(result)) {
        cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
        MLE_SE = result
      } else {
        MLE_SE <- result
      }
      MLE <- data.table("param" = names(par0_constrained),
                        "estimates" = MLE_par,
                        "sd(+/-)" = MLE_SE)
      
      print("The retrieved parameters are:")
      prmatrix(MLE)
      
      if (MLE_estimates$convergence == 0) {
        convergence <- "TRUE"
        cat(paste0("\033[0;32m","#################### CONVERGENCE: GLOBAL #####################","\033[0m","\n"))
      } else {
        convergence = "FALSE"
        cat(paste0("\033[0;34m","#################### CONVERGENCE: LOCAL #####################","\033[0m","\n"))
      }
      end.time = Sys.time(); time_taken <- end.time - start.time
      return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
      
      cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
      
  }  
  
  #==========================================================================
  # For Unconstrained full inversion
  #==========================================================================
  if (unconstrained == TRUE) {
    
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    Sys.sleep(1)
    if (obj.fn == "log-LL") {
      print("Log-Likelihood will be used to construct objective function")
      ##Create log-likelihood function
      
      NLL_unconstr = function(pars, data) {
        
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
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 = pars[3],
            
            z = pars[4],
            rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        # Negative log-likelihood
        smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                           log = TRUE))
        return(smull)
      }
    } else {
      
      print("Spectral Error index from Lee et al. 1999 will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_unconstr = function(pars, data) {
        
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
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            #a_non_water_path = IOP_files[idx_a],
            #bb_non_water_path = IOP_files[idx_bb],
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 = pars[3],
            
            z = pars[4],
            rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
            
            
            Rrs_input_for_slope = data,
            
            slope.parametric = auto_spectral_slope,
            
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        rrs_est = Gpred[[1]]$Rrs
        
        # Sum-squared residual of error (SSE)
        #sse= sum((data - rrs_est)^2)
        #return(sse)
        
        #The Spectral error index from Lee 1999
        # Define the spectral regions
        region1 <- which(wavelength >= 400 & wavelength <= 675)
        region2 <- which(wavelength >= 750 & wavelength <= 830)
        
        # Calculate the numerator of the error index
        numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                            sum((data[region2] - rrs_est[region2])^2))
        
        # Calculate the denominator of the error index
        denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
        
        # Calculate the error index
        smull <- numerator / denominator
        
        return(smull)
        
      }
      
    }
    
    ##Optimize the error function
    cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
    
    if (auto_spectral_slope == TRUE) {
      print("Spectral slopes will be calculated using QAA")
    } else {
      if (manual_spectral_slope == TRUE) {
        print(paste0("Spectral slopes will be used from user supplied values: ",
                     manual_spectral_slope_vals[1], ", " , manual_spectral_slope_vals[2], ", " ,manual_spectral_slope_vals[3]))
      } else {
        print("Spectral slopes will be used same as model default")
      }
    }
    
    #Instantiate initial values
    par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], 
             
             z = initial[4], 
             initial[5:(4+initial_rb_length)],
             pop.sd = initial[length(initial)])
    
    
    cat(paste0("\033[0;33m","####################UNCONSTRAINED INVERSION#########################","\033[0m","\n"))
    
    cat(paste0("\033[0;32m","Initial values are: [chl]=",par0[1], ", adg(440)=",par0[2], ", bbp(555)=", par0[3], ", zB=", par0[4],",  RB={", toString(as.numeric(par0[5:(length(par0)-1)])),"},  population.sigma=", par0[length(par0)],"\033[0m","\n"))
    
    
    cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
    #Sys.sleep(2)
    start.time = Sys.time()
    
    if (sa.model == "am03") {
      print("Albert & Mobley 2003 SA model used for forward model")
    } else {
      print("Lee et al. 1999 SA model used for forward model")
    }
    
    if (method.opt == "L-BFGS-B") {
      
      MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata, 
                            lower = lower.b,     # Lower bound on parameters
                            upper = upper.b,  # Upper bound on parameters
                            #gr = grad.calculate,
                            method = method.opt,
                            control = list(#fnscale = 1, 
                              parscale = abs(par0)),
                            #hessian = TRUE
      )
    } 
    if (method.opt == "levenberg-marqardt") {
      
      LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_unconstr,data = obsdata, print.info = F)
      MLE_estimates = data.frame("par"=LM_estimates$b)
      
    } 
    
    if (method.opt == "auglag") {
      print("Augmented Lagriangian with equality constraints will be used for inversion")
      
      
      fheq <- function(pars,data) sum(pars[5:(length(pars)-1)]) - 1
      
      fhin <- function(pars,data) c(
        #bounds for aerial fractions
        pars[5:(length(pars)-1)]
      )
      
      MLE_estimates <- alabama::auglag(fn = NLL_unconstr,par = par0, 
                                       #gr = pracma::grad(NLL_unconstr, x0 = par, data= obsdata),
                                       heq = fheq,
                                       hin = fhin,
                                       #lower = lower.bound,
                                       #upper = upper.bound,
                                       #data=surface_rrs_translate(insitu.data),
                                       data = obsdata,
                                       control.outer = list(trace = F, method = "nlminb")
      )
      
      #print(inverse_output_fmincon$par, digits=3)
      
      #print(sum(inverse_output_fmincon$par[5:9]))
      
      #sprintf("%.5f", inverse_output_fmincon$par)
      
    }
    
    
    if (method.opt == "Nelder-Mead" |method.opt ==  "SANN" | method.opt ==  "Brent") {
      
      print("Simplex based Optimization will be used for inversion")
      MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata, 
                            #lower = c(0, 0, 0),     # Lower bound on parameters
                            #upper = c(10, 2, 0.3),  # Upper bound on parameters
                            method = method.opt,
                            control = list(#fnscale = 1, 
                              parscale = as.numeric(abs(par0))),
                            hessian = FALSE)
      
      
      
      
    }
    
    
    cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
    
    #print(MLE_estimates$par)
    
    cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
    
    #Calculate hessian matrix for var-covar matrix
    if (method.opt == "auglag") {
      hessian.inverse <- MLE_estimates$hessian 
    } else {
      hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_unconstr, data=obsdata)
    }
    #hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_unconstr, data=obsdata)
    rownames(hessian.inverse) <- names(par0)
    colnames(hessian.inverse) <- names(par0)
    cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
    prmatrix(hessian.inverse)
    cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
    
    MLE_par <- MLE_estimates$par
    
    result <- tryCatch({
      sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
      error = myFun
    )
    if (!is.numeric(result)) {
      cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
      MLE_SE = result
    } else {
      MLE_SE <- result
    }
    MLE <- data.table("param" = names(par0),
                      "estimates" = MLE_par,
                      "sd(+/-)" = MLE_SE)
    
    print("The retrieved parameters are:")
    prmatrix(MLE)
    
    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
      cat(paste0("\033[0;32m","#################### CONVERGENCE: GLOBAL #####################","\033[0m","\n"))
    } else {
      convergence = "FALSE"
      cat(paste0("\033[0;34m","#################### CONVERGENCE: LOCAL #####################","\033[0m","\n"))
    }
    end.time = Sys.time(); time_taken <- end.time - start.time
    return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
    
    cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
    
    
    
  }
  
  }
  
#===============================================================================================================
#Solver for deep waters
#===============================================================================================================
solve.objective.inverse.deep.final.fast <- function(
    bbp.constrain = F,
    bbp.constrain.value,
    
    auto_spectral_slope = T,
    manual_spectral_slope = F,
    manual_spectral_slope_vals = c("s_g"=0.015, "s_d"=0.01160, "gamma"=0.5),
    
    obj.fn, initial, sa.model="am03", method.opt,
    
    obsdata,
    lower.b, upper.b, wave = wavelength 
    #,batch, pop.sd
    ){
  if (auto_spectral_slope == TRUE & manual_spectral_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  
  myFun <- function(x) {
    NA
  }
  #initial <- abs(initial)
  
  if (bbp.constrain == TRUE) {
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### bbp CONSTRAINED deep WATER INVERSION #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    
    if (obj.fn == "log-LL") {
      
      print("Log-Likelihood will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =bbp.constrain.value,
            
            z = zB,
            rb.fraction = fA.set,
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =bbp.constrain.value,
            
            z = zB,
            rb.fraction = fA.set,
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        # Negative log-likelihood
        
        smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], log = TRUE)) 
        return(smull)
      }
    } 
    if (obj.fn == "obj_L98") {
      print("Spectral error index from Lee et al. 1999 will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =bbp.constrain.value,
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =bbp.constrain.value,
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        rrs_est = Gpred[[1]]$Rrs
        
        # Sum-squared residual of error (SSE)
        #sse= sum((data - rrs_est)^2)
        #return(sse)
        
        #The Spectral error index from Lee 1999
        # Define the spectral regions
        region1 <- which(wavelength >= 400 & wavelength <= 675)
        region2 <- which(wavelength >= 750 & wavelength <= 830)
        
        # Calculate the numerator of the error index
        numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                            sum((data[region2] - rrs_est[region2])^2))
        
        # Calculate the denominator of the error index
        denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
        
        # Calculate the error index
        err <- numerator / denominator
        
        return(err)
      }
    }
    
    ##Optimize the log-likelihood function
    cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
    
    if (auto_spectral_slope == TRUE) {
      print("Spectral slopes will be calculated using QAA")
    } else {
      if (manual_spectral_slope == TRUE) {
        print(paste0("Spectral slopes will be used from user supplied values: ",
                     manual_spectral_slope_vals[1], ", " , manual_spectral_slope_vals[2], ", " ,manual_spectral_slope_vals[3]))
      } else {
        print("Spectral slopes will be used same as model default")
      }
    }
    
    par0 = c(chl = initial[1], adg440 = initial[2], pop.sd = initial[length(initial)])
    cat(paste0("\033[0;36m","Initial values are ",par0[1],"    ", par0[2],"\033[0m","\n"))
    
    
    cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
    #Sys.sleep(2)
    start.time = Sys.time()
    
    if (sa.model == "am03") {
      print("Albert & Mobley 2003 SA model used for forward model")
    } else {
      print("Lee et al. 1999 SA model used for forward model")
    }
    
    if (method.opt == "L-BFGS-B") {
      #browser()
      MLE_estimates = optim(par = par0, fn = NLL_deep, data = obsdata, 
                            lower = lower.b,     # Lower bound on parameters
                            upper = upper.b,  # Upper bound on parameters
                            #gr = saber.forward.jacobian.analytical.optim,
                            # chl = par[1],
                            # acdom.440 = par[2],
                            # anap.440 = par[3],
                            method = method.opt,
                            control = list(parscale = abs(par0)),
                            hessian = F)
    } else {
      if (method.opt == "levenberg-marqardt") {
        
        LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_deep,data = obsdata, print.info = F)
        MLE_estimates = data.frame("par"=LM_estimates$b)
        
      } else {
        
        MLE_estimates = optim(par = par0, fn = NLL_deep, data = obsdata, 
                              #lower = c(0, 0, 0),     # Lower bound on parameters
                              #upper = c(10, 2, 0.3),  # Upper bound on parameters
                              method = method.opt,
                              control = list(parscale = abs(par0)),
                              hessian = FALSE)
      }
      
    }
    
    cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
    
    #print(MLE_estimates$par)
    
    cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
    #Sys.sleep(2)
    
    #Calculate hessian matrix for var-covar matrix
    hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_deep, data=obsdata)
    rownames(hessian.inverse) <- names(par0)
    colnames(hessian.inverse) <- names(par0)
    cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
    prmatrix(hessian.inverse)
    cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
    
    MLE_par <- MLE_estimates$par
    result <- tryCatch({
      sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
      error = myFun
    )
    if (!is.numeric(result)) {
      cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
      MLE_SE = result
    } else {
      MLE_SE <- result
    }
    #MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
    MLE <- data.table("param" = names(par0),
                      "estimates" = MLE_par,
                      "sd(+/-)" = MLE_SE)
    
    print("The retrieved parameters are:")
    prmatrix(MLE)
    cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
    } else {
      convergence = "FALSE"
    }
    end.time = Sys.time(); time_taken <- end.time - start.time
    return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
   
  } else {
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### UNCONSTRAINED DEEP WATER INVERSION #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    if (obj.fn == "log-LL") {
      
      print("Log-Likelihood will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =pars[3],
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =pars[3],
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        # Negative log-likelihood
        
        smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], log = TRUE)) 
        return(smull)
      }
    } 
    if (obj.fn == "obj_L98") {
      print("Spectral error index from Lee et al. 1999 will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =pars[3],
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            a_non_water_path = abs_path,
            bb_non_water_path = bb_path,
            
            chl = pars[1], 
            a_dg = pars[2],
            bbp.550 =pars[3],
            
            z = zB,
            rb.fraction = fA.set,
            
            
            slope.parametric = auto_spectral_slope,
            Rrs_input_for_slope = data, 
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            verbose = F, wavelength = wave
          )
        }
        
        rrs_est = Gpred[[1]]$Rrs
        
        # Sum-squared residual of error (SSE)
        #sse= sum((data - rrs_est)^2)
        #return(sse)
        
        #The Spectral error index from Lee 1999
        # Define the spectral regions
        region1 <- which(wavelength >= 400 & wavelength <= 675)
        region2 <- which(wavelength >= 750 & wavelength <= 830)
        
        # Calculate the numerator of the error index
        numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                            sum((data[region2] - rrs_est[region2])^2))
        
        # Calculate the denominator of the error index
        denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
        
        # Calculate the error index
        err <- numerator / denominator
        
        return(err)
      }
    }
    
    ##Optimize the log-likelihood function
    cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
    
    if (auto_spectral_slope == TRUE) {
      print("Spectral slopes will be calculated using QAA")
    } else {
      if (manual_spectral_slope == TRUE) {
        print(paste0("Spectral slopes will be used from user supplied values: ",
                     manual_spectral_slope_vals[1], ", " , manual_spectral_slope_vals[2], ", " ,manual_spectral_slope_vals[3]))
      } else {
        print("Spectral slopes will be used same as model default")
      }
    }
    
    par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], pop.sd = initial[length(initial)])
    cat(paste0("\033[0;36m","Initial values are [chl]=",par0[1],"   [adg440]= ", par0[2],"   [bbp555]= ", par0[3],"\033[0m","\n"))
    
    
    cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
    #Sys.sleep(2)
    start.time = Sys.time()
    
    if (sa.model == "am03") {
      print("Albert & Mobley 2003 SA model used for forward model")
    } else {
      print("Lee et al. 1999 SA model used for forward model")
    }
    
    if (method.opt == "L-BFGS-B") {
      #browser()
      MLE_estimates = optim(par = par0, fn = NLL_deep, data = obsdata, 
                            lower = lower.b,     # Lower bound on parameters
                            upper = upper.b,  # Upper bound on parameters
                            #gr = saber.forward.jacobian.analytical.optim,
                            # chl = par[1],
                            # acdom.440 = par[2],
                            # anap.440 = par[3],
                            method = method.opt,
                            control = list(parscale = abs(par0)),
                            hessian = F)
    } else {
      if (method.opt == "levenberg-marqardt") {
        
        LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_deep,data = obsdata, print.info = F)
        MLE_estimates = data.frame("par"=LM_estimates$b)
        
      } else {
        
        MLE_estimates = optim(par = par0, fn = NLL_deep, data = obsdata, 
                              #lower = c(0, 0, 0),     # Lower bound on parameters
                              #upper = c(10, 2, 0.3),  # Upper bound on parameters
                              method = method.opt,
                              control = list(parscale = abs(par0)),
                              hessian = FALSE)
      }
      
    }
    
    cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
    
    #print(MLE_estimates$par)
    
    cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
    #Sys.sleep(2)
    
    #Calculate hessian matrix for var-covar matrix
    hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_deep, data=obsdata)
    rownames(hessian.inverse) <- names(par0)
    colnames(hessian.inverse) <- names(par0)
    cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
    prmatrix(hessian.inverse)
    cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
    
    MLE_par <- MLE_estimates$par
    result <- tryCatch({
      sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
      error = myFun
    )
    if (!is.numeric(result)) {
      cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
      MLE_SE = result
    } else {
      MLE_SE <- result
    }
    #MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
    MLE <- data.table("param" = names(par0),
                      "estimates" = MLE_par,
                      "sd(+/-)" = MLE_SE)
    
    print("The retrieved parameters are:")
    prmatrix(MLE)
    cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
    } else {
      convergence = "FALSE"
    }
    end.time = Sys.time(); time_taken <- end.time - start.time
    return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
    
  }
  
  
}
