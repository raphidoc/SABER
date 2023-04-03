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

solve.objective.inverse.shallow.final <- function(
                                          constrained = FALSE,
                                          constrain.param, 
                                          bbp.550 =0.05,
                                          
                                          auto_spectral_slope = T,
                                          manual_spectral_slope = F,
                                          manual_spectral_slope_vals = c("s_g"=0.015, "s_d"=0.01160, "gamma"=0.5),
                                          
                                          obj.fn, initial, sa.model="am03", method.opt,
                                          
                                          obsdata,  
                                          
                                          lower.b, upper.b, 
                                          batch=FALSE, pop.sd=FALSE){
  
  if (auto_spectral_slope == TRUE & manual_spectral_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  myFun <- function(x) {
    NA
  }
  
  if (constrained == TRUE) {
    
    #initial <- abs(initial)
    if (obj.fn == "log-LL") {
      print("Log-Likelihood will be used to construct objective function")
      
      ##Create log-likelihood function
      NLL_constr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_final(
                                use_true_IOPs = F, 
                                a_non_water_path = IOP_files[idx_a],
                                bb_non_water_path = IOP_files[idx_bb],
                                
                                chl = constrain.param[1], 
                                acdom440 =NULL, 
                                anap440 =NULL , 
                                a_dg = constrain.param[2],
                                bbp.550 = bbp.550,
                                
                                z = pars[1],
                                use_spectral_rb = F, 
                                rb.fraction = as.numeric(pars[2:6]),
                                
                  
                                realdata = data,
                                
                                slope.parametric = auto_spectral_slope,
                                dg_composite = T,
                                
                                use_manual_slope =manual_spectral_slope,
                                manual_slope =  manual_spectral_slope_vals,
                                
                                use_spectral_shape_chl = F,
                                use_spectral_shape_dg = T,
                                
                                sicf = F, q_phi = 0.02, 
                                fDOM = F,
                                verbose = F, plot = F
                                )
          
          
        } else {
          Gpred = Lee_forward(chl = constrain.param[1], acdom440 = constrain.param[2],
                              anap440 = constrain.param[3], bbp.550 = bbp.550,
                              z = pars[1], rb.fraction = as.numeric(pars[2:6]),
                              verbose = F, realdata = data, plot = F)
        }
        
        # Negative log-likelihood
        if(pop.sd == "TRUE") {
          if (batch == FALSE) {
            
            #ON for Albert Rrs
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) 
            
          } else {
            
            #ON for IOCCG
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, 
                               log = TRUE)) 
            
          }
        } else {
          if (batch == TRUE | batch == FALSE) {
            
            #ON for unknown sigma
            smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[7], 
                               log = TRUE)) 
            
          }
          
        }
        return(smull)
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
      
      #Turn ON known population sigma for likelihood
      if(pop.sd == TRUE) {
        par0 = c(chl = initial[1], adg440 = initial[2], #anap440 = initial[3],
                 z = initial[3], 
                 #rb0 = initial[5], 
                 rb1 = initial[4], rb2 = initial[5],
                 rb3 = initial[6], rb4 = initial[7], 
                 rb5 = initial[8]
        )#, bbp550 = initial[4])
      } else {
        #Turn ON for unknown population sigma for likelihood
        par0 = c(chl = initial[1], adg440 = initial[2], #anap440 = initial[3], 
                 z = initial[3], 
                 #rb0 = initial[5], 
                 rb1 = initial[4], rb2 = initial[5],
                 rb3 = initial[6], rb4 = initial[7], 
                 rb5 = initial[8],
                 pop.sd = initial[9])
      }
      par0_constrained = par0[-(1:2)]
      
      
      cat(paste0("\033[0;33m","####################CONSTRAINED INVERSION#########################","\033[0m","\n"))
      
      cat(paste0("\033[0;32m","Initial values are: zB=", par0[4],"  RB[1:5]={", par0[5],",", par0[6],",", par0[7],",", par0[8],",", par0[9],"}  population.sigma=", par0[10],"\033[0m","\n"))
      
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      Sys.sleep(2)
      start.time = Sys.time()
      
      if (sa.model == "am03") {
        print("Albert & Mobley 2003 SA model used for forward model")
      } else {
        print("Lee et al. 1999 SA model used for forward model")
      }
      
      if (method.opt == "L-BFGS-B") {
        
        print("L-BFGS-B Optimization will be used for inversion")
        MLE_estimates = optim(par = as.numeric(par0_constrained), fn = NLL_constr, data = obsdata, 
                              lower = as.numeric(lower.b[-(1:2)]),     # Lower bound on parameters
                              upper = as.numeric(upper.b[-{1:2}]),  # Upper bound on parameters
                              #gr = grad.calculate,
                              method = method.opt,
                              control = list(#fnscale = 1, 
                                parscale = as.numeric(abs(par0_constrained))),
                              #hessian = TRUE
        )
      } else {
        if (method.opt == "levenberg-marqardt") {
          
          print("Levenberg-Marquardt Optimization will be used for inversion")
          LM = marqLevAlg::marqLevAlg(b = par0_constrained, fn = NLL_constr,data = obsdata, print.info = F)
          MLE_estimates = data.frame("par"=LM_estimates$b)
          
        } else {
          print("Simplex based Optimization will be used for inversion")
          MLE_estimates = optim(par = par0_constrained, fn = NLL_constr, data = obsdata, 
                                #lower = c(0, 0, 0),     # Lower bound on parameters
                                #upper = c(10, 2, 0.3),  # Upper bound on parameters
                                method = method.opt,
                                control = list(parscale = abs(par0)),
                                hessian = FALSE)
        }
        
      }
      
      if (method.opt == "auglag") {
        print("Augmented Lagriangian with equality constraints will be used for inversion")
        
        if (pop.sd == TRUE) {
          fheq <- function(pars,data) sum(pars[2:(length(pars))]) - 1
          
          fhin <- function(pars,data) c(#bounds for [chl]
            # 0.5 - pars[1],
            # pars[1] - 10,
            # 
            # #bounds for adg443
            # 0.1 - pars[2],
            # pars[2] - 10,
            # 
            # #bounds for bbp555
            # 0.0001 - pars[3],
            # pars[3] - 0.05,
            # 
            # #bounds for zB
            # 0.1 - pars[4],
            # pars[4] - 10,
            
            #bounds for aerial fractions
            pars[2:(length(pars))]
            #,
            
            # #bounds for zB
            # 0.0001 - pars[10],
            # pars[10] - 1
            
          )
          
        } else {
          fheq <- function(pars,data) sum(pars[2:(length(pars)-1)]) - 1
          
          fhin <- function(pars,data) c(
            #bounds for [chl]
            # 0.5 - pars[1],
            # pars[1] - 10,
            # 
            # #bounds for adg443
            # 0.1 - pars[2],
            # pars[2] - 10,
            # 
            # #bounds for bbp555
            # 0.0001 - pars[3],
            # pars[3] - 0.05,
            # 
            # #bounds for zB
            # 0.1 - pars[4],
            # pars[4] - 10,
            
            #bounds for aerial fractions
            pars[2:(length(pars)-1)]
            #,
            
            # #bounds for pop sd
            # 0.0001 - pars[10],
            # pars[10] - 1
            
          )
        }
        
        
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
      
    } else {
      #@@@@@@@@@@@@@@ SSR Objective function @@@@@@@@@@@@@@@@@@@@@@@@
      if (obj.fn == "obj_L98") {
        print("Spectral Error index from Lee et al. 1999 will be used to construct objective function")
        
        ##Create log-likelihood function
        NLL_constr = function(pars, data) {
          
          # Values predicted by the forward model for single RUN
          if (sa.model == "am03") {
                                Gpred = Saber_forward_final(
                                  use_true_IOPs = F, 
                                  a_non_water_path = IOP_files[idx_a],
                                  bb_non_water_path = IOP_files[idx_bb],
                                  
                                  chl = constrain.param[1], 
                                  acdom440 =NULL, 
                                  anap440 =NULL , 
                                  a_dg = constrain.param[2],
                                  bbp.550 = bbp.550,
                                  
                                  z = pars[1],
                                  use_spectral_rb = F, 
                                  rb.fraction = as.numeric(pars[2:6]),
                                  
                                  
                                  realdata = data,
                                  
                                  slope.parametric = auto_spectral_slope,
                                  dg_composite = T,
                                  
                                  use_manual_slope =manual_spectral_slope,
                                  manual_slope =  manual_spectral_slope_vals,
                                  
                                  use_spectral_shape_chl = F,
                                  use_spectral_shape_dg = T,
                                  
                                  sicf = F, q_phi = 0.02, 
                                  fDOM = F,
                                  verbose = F, plot = F
            )
            
            
          } else {
            Gpred = Lee_forward(chl = constrain.param[1], acdom440 = constrain.param[2],
                                anap440 = constrain.param[3], bbp.550 = bbp.550,
                                z = pars[1], rb.fraction = as.numeric(pars[2:6]),
                                verbose = F, realdata = data, plot = F)
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
        
        #Turn ON known population sigma for likelihood
        if(pop.sd == TRUE) {
          par0 = c(chl = initial[1], adg440 = initial[2], #anap440 = initial[3],
                   z = initial[3], 
                   #rb0 = initial[5], 
                   rb1 = initial[4], rb2 = initial[5],
                   rb3 = initial[6], rb4 = initial[7], 
                   rb5 = initial[8]
          )#, bbp550 = initial[4])
        } else {
          #Turn ON for unknown population sigma for likelihood
          par0 = c(chl = initial[1], adg440 = initial[2], #anap440 = initial[3], 
                   z = initial[3], 
                   #rb0 = initial[5], 
                   rb1 = initial[4], rb2 = initial[5],
                   rb3 = initial[6], rb4 = initial[7], 
                   rb5 = initial[8],
                   pop.sd = initial[9])
        }
        par0_constrained = par0[-(1:2)]
        
        
        cat(paste0("\033[0;33m","####################CONSTRAINED INVERSION#########################","\033[0m","\n"))
        
        cat(paste0("\033[0;32m","Initial values are: zB=", par0[4],"  RB[1:5]={", par0[5],",", par0[6],",", par0[7],",", par0[8],",", par0[9],"}  population.sigma=", par0[10],"\033[0m","\n"))
        
        cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
        Sys.sleep(2)
        start.time = Sys.time()
        
        if (sa.model == "am03") {
          print("Albert & Mobley 2003 SA model used for forward model")
        } else {
          print("Lee et al. 1999 SA model used for forward model")
        }
        
        if (method.opt == "L-BFGS-B") {
          
          MLE_estimates = optim(par = as.numeric(par0_constrained), fn = NLL_constr, data = obsdata, 
                                lower = as.numeric(lower.b[-(1:2)]),     # Lower bound on parameters
                                upper = as.numeric(upper.b[-{1:2}]),  # Upper bound on parameters
                                #gr = grad.calculate,
                                method = method.opt,
                                control = list(#fnscale = 1, 
                                  parscale = as.numeric(abs(par0_constrained))),
                                #hessian = TRUE
          )
        } else {
          if (method.opt == "levenberg-marqardt") {
            
            LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_constr,data = obsdata, print.info = F, minimize = T)
            MLE_estimates = data.frame("par"=LM_estimates$b)
            
          } else {
            
            MLE_estimates = optim(par = par0, fn = NLL_constr, data = obsdata, 
                                  #lower = c(0, 0, 0),     # Lower bound on parameters
                                  #upper = c(10, 2, 0.3),  # Upper bound on parameters
                                  method = method.opt,
                                  control = list(parscale = abs(par0)),
                                  hessian = FALSE)
          }
          
        }
        
        if (method.opt == "auglag") {
          print("Augmented Lagriangian with equality constraints will be used for inversion")
          
          if (pop.sd == TRUE) {
            fheq <- function(pars,data) sum(pars[2:(length(pars))]) - 1
            
            fhin <- function(pars,data) c(#bounds for [chl]
              # 0.5 - pars[1],
              # pars[1] - 10,
              # 
              # #bounds for adg443
              # 0.1 - pars[2],
              # pars[2] - 10,
              # 
              # #bounds for bbp555
              # 0.0001 - pars[3],
              # pars[3] - 0.05,
              # 
              # #bounds for zB
              # 0.1 - pars[4],
              # pars[4] - 10,
              
              #bounds for aerial fractions
              pars[2:(length(pars))]
              #,
              
              # #bounds for zB
              # 0.0001 - pars[10],
              # pars[10] - 1
              
            )
            
          } else {
            fheq <- function(pars,data) sum(pars[2:(length(pars)-1)]) - 1
            
            fhin <- function(pars,data) c(#bounds for [chl]
              # 0.5 - pars[1],
              # pars[1] - 10,
              # 
              # #bounds for adg443
              # 0.1 - pars[2],
              # pars[2] - 10,
              # 
              # #bounds for bbp555
              # 0.0001 - pars[3],
              # pars[3] - 0.05,
              # 
              # #bounds for zB
              # 0.1 - pars[4],
              # pars[4] - 10,
              
              #bounds for aerial fractions
              pars[2:(length(pars)-1)]
              #,
              
              # #bounds for zB
              # 0.0001 - pars[10],
              # pars[10] - 1
              
            )
          }
          
          
          MLE_estimates <- alabama::auglag(fn = NLL_constr,par = par0[-length(par0)], 
                                           #gr = pracma::grad(NLL_unconstr, x0 = par, data= obsdata),
                                           heq = fheq,
                                           hin = fhin,
                                           #lower = lower.bound,
                                           #upper = upper.bound,
                                           #data=surface_rrs_translate(insitu.data),
                                           data = obsdata,
                                           control.outer = list(trace = F, method = "nlminb")
          )
          
          
          #print(MLE_estimates$par, digits=3)
          
          #print(sum(MLE_estimates$par[5:9]))
          
          #sprintf("%.5f", MLE_estimates$par)
          
        }
        cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
        
        #print(MLE_estimates$par)
        
        cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
        Sys.sleep(2)
        
        #Calculate hessian matrix for var-covar matrix
        if (method.opt == "auglag") {
          hessian.inverse <- MLE_estimates$hessian 
        } else {
          hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_constr, data=obsdata)
        }
        #hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL_constr, data=obsdata)
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
      
      if (obj.fn == "SSR") {
        print("Sum-square of residual will be used to construct objective function")
        
        cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
        Fit <- data.frame("C_ph"=initial[1],
                          "a_cdom.440"=initial[2],
                          "a.nap.440"=initial[3])
        
        params <- as.numeric(Fit)
        cat(paste0("\033[0;36m","Initial values are ",params[1],"    ", params[2],"    ", 
                   params[3],"\033[0m","\n"))
        cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
        Sys.sleep(1)
        start.time = Sys.time()
        
        ##Single RUN
        inverse_output <- pracma::fminunc(fn = Saber_inverse,x0 = params,valdata=obsdata,
                                          bbp.550= Fit.input$bbp.550)
        
        cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
        
        Fit.optimized.ssobj <- data.frame("chl"=inverse_output$par[1], 
                                          "acdom.440"=inverse_output$par[2],
                                          "anap.440"=inverse_output$par[3])
        prmatrix(Fit.optimized.ssobj)
        
        #Calculate hessian matrix for var-covar matrix
        hessian.inverse <- numDeriv::hessian(x =as.numeric(Fit.optimized.ssobj), 
                                             func = Saber_inverse, valdata=obsdata,
                                             bbp.550= Fit.input$bbp.550)
        rownames(hessian.inverse) <- names(par0)
        colnames(hessian.inverse) <- names(par0)
        cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
        prmatrix(hessian.inverse)
        cat(paste0("\033[0;36m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
        
        MLE_par <- as.numeric(Fit.optimized.ssobj)
        MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
        MLE <- data.table("param" = names(par0),
                          "estimates" = MLE_par,
                          "sd(+/-)" = MLE_SE)
        
        print("The retrieved parameters are:")
        prmatrix(MLE)
        cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
        end.time = Sys.time(); time_taken <- end.time - start.time
        return(list(Fit.optimized.ssobj,"time.elapsed"= time_taken))
      }
    }
    
  } else {
    cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
    cat(paste0("\033[0;39m","########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######","\033[0m","\n"))
    cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
    Sys.sleep(2)
    if (obj.fn == "log-LL") {
      print("Log-Likelihood will be used to construct objective function")
      ##Create log-likelihood function
      
      NLL_unconstr = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_final(
                                use_true_IOPs = F, 
                                a_non_water_path = IOP_files[idx_a],
                                bb_non_water_path = IOP_files[idx_bb],
                                
                                chl = pars[1], 
                                acdom440 =NULL, 
                                anap440 =NULL , 
                                a_dg = pars[2],
                                bbp.550 = pars[3],
                                
                                z = pars[4],
                                use_spectral_rb = F, 
                                rb.fraction = pars[5:(length(pars)-1)],
                                
                  
                                realdata = data,
                                
                                slope.parametric = auto_spectral_slope,
                                dg_composite = T,
                                
                                use_manual_slope =manual_spectral_slope,
                                manual_slope =  manual_spectral_slope_vals,
                                
                                use_spectral_shape_chl = F,
                                use_spectral_shape_dg = T,
                                
                                sicf = F, q_phi = 0.02, 
                                fDOM = F,
                                verbose = F, plot = F
                                )
          
        } else {
          Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                              anap440 =pars[3], bbp.550 = bbp.550,
                              z = pars[4], rb.fraction = pars[5:(length(pars)-1)],
                              verbose = F, realdata = data, plot = F)
        }
        
        # Negative log-likelihood
        if(pop.sd == "TRUE") {
          if (batch == FALSE) {
            
            #ON for Albert Rrs
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) 
            
          } else {
            
            #ON for IOCCG
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, 
                               log = TRUE)) 
            
          }
        } else {
          if (batch == TRUE | batch == FALSE) {
            
            #ON for unknown sigma
            smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                               log = TRUE)) 
            
          }
          
        }
        return(smull)
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
      #Turn ON known population sigma for likelihood
      if(pop.sd == TRUE) {
        par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3],
                 z = initial[4], 
                 #rb0 = initial[5], 
                 rb1 = initial[5], rb2 = initial[6],
                 rb3 = initial[7], rb4 = initial[8], 
                 rb5 = initial[9]
        )#, bbp550 = initial[4])
      } else {
        #Turn ON for unknown population sigma for likelihood
        par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], 
                 z = initial[4], 
                 #rb0 = initial[5], 
                 rb1 = initial[5], rb2 = initial[6],
                 rb3 = initial[7], rb4 = initial[8], 
                 rb5 = initial[9],
                 pop.sd = initial[10])
      }
      
      cat(paste0("\033[0;36m","Initial values are: chl=",par0[1],", adg[440]=", par0[2], ", bbp[550]=", par0[3], 
                  ", zB=", par0[4],", RB[1:5]={", par0[5]," ", par0[6]," ", par0[7]," ", par0[8]," ", par0[9],"}, population.sigma=", par0[10],"\033[0m","\n"))
      
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      
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
        
        if (pop.sd == TRUE) {
          fheq <- function(pars,data) sum(par0[5:(length(pars))]) - 1
          
          fhin <- function(pars,data) c(#bounds for [chl]
            # 0.5 - pars[1],
            # pars[1] - 10,
            # 
            # #bounds for adg443
            # 0.1 - pars[2],
            # pars[2] - 10,
            # 
            # #bounds for bbp555
            # 0.0001 - pars[3],
            # pars[3] - 0.05,
            # 
            # #bounds for zB
            # 0.1 - pars[4],
            # pars[4] - 10,
            
            #bounds for aerial fractions
            pars[5:(length(pars))]
            #,
            
            # #bounds for zB
            # 0.0001 - pars[10],
            # pars[10] - 1
            
          )
          
        } else {
          fheq <- function(pars,data) sum(pars[5:(length(pars)-1)]) - 1
          
          fhin <- function(pars,data) c(#bounds for [chl]
            # 0.5 - pars[1],
            # pars[1] - 10,
            # 
            # #bounds for adg443
            # 0.1 - pars[2],
            # pars[2] - 10,
            # 
            # #bounds for bbp555
            # 0.0001 - pars[3],
            # pars[3] - 0.05,
            # 
            # #bounds for zB
            # 0.1 - pars[4],
            # pars[4] - 10,
            
            #bounds for aerial fractions
            pars[5:(length(pars)-1)]
            #,
            
            # #bounds for zB
            # 0.0001 - pars[10],
            # pars[10] - 1
            
          )
        }
        
        
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
        
        
        #print(MLE_estimates$par, digits=3)
        
        #print(sum(MLE_estimates$par[5:9]))
        
        #sprintf("%.5f", MLE_estimates$par)
        
      } else {
        
        MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata, 
                              #lower = c(0, 0, 0),     # Lower bound on parameters
                              #upper = c(10, 2, 0.3),  # Upper bound on parameters
                              method = method.opt,
                              control = list(parscale = abs(par0)),
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
      
    } else {
      if (obj.fn == "obj_L98") {
        
        print("Spectral error index from Lee et al. 1999 will be used to construct objective function")
        
        ##Create log-likelihood function
        
        NLL_unconstr = function(pars, data) {
          
          # Values predicted by the forward model for single RUN
          if (sa.model == "am03") {
            Gpred = Saber_forward_final(
              use_true_IOPs = F, 
              a_non_water_path = IOP_files[idx_a],
              bb_non_water_path = IOP_files[idx_bb],
              
              chl = pars[1], 
              acdom440 =NULL, 
              anap440 =NULL , 
              a_dg = pars[2],
              bbp.550 = pars[3],
              
              z = pars[4],
              use_spectral_rb = F, 
              rb.fraction = pars[5:(length(pars)-1)],
              
              
              realdata = data,
              
              slope.parametric = auto_spectral_slope,
              dg_composite = T,
              
              use_manual_slope =manual_spectral_slope,
              manual_slope =  manual_spectral_slope_vals,
              
              use_spectral_shape_chl = F,
              use_spectral_shape_dg = T,
              
              sicf = F, q_phi = 0.02, 
              fDOM = F,
              verbose = F, plot = F
            )
            
          } else {
            Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                                anap440 =pars[3], bbp.550 = bbp.550,
                                z = pars[4], rb.fraction = pars[5:(length(pars)-1)],
                                verbose = F, realdata = data, plot = F)
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
        #Turn ON known population sigma for likelihood
        if(pop.sd == TRUE) {
          par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3],
                   z = initial[4], 
                   #rb0 = initial[5], 
                   rb1 = initial[5], rb2 = initial[6],
                   rb3 = initial[7], rb4 = initial[8], 
                   rb5 = initial[9]
          )#, bbp550 = initial[4])
        } else {
          #Turn ON for unknown population sigma for likelihood
          par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], 
                   z = initial[4], 
                   #rb0 = initial[5], 
                   rb1 = initial[5], rb2 = initial[6],
                   rb3 = initial[7], rb4 = initial[8], 
                   rb5 = initial[9],
                   pop.sd = initial[10])
        }
        
        cat(paste0("\033[0;36m","Initial values are: chl=",par0[1],", adg[440]=", par0[2], ", bbp[550]=", par0[3], 
                   ", zB=", par0[4],", RB[1:5]={", par0[5]," ", par0[6]," ", par0[7]," ", par0[8]," ", par0[9],"}, population.sigma=", par0[10],"\033[0m","\n"))
        
        cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
        
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
          
          if (pop.sd == TRUE) {
            fheq <- function(pars,data) sum(par0[5:(length(pars))]) - 1
            
            fhin <- function(pars,data) c(#bounds for [chl]
              # 0.5 - pars[1],
              # pars[1] - 10,
              # 
              # #bounds for adg443
              # 0.1 - pars[2],
              # pars[2] - 10,
              # 
              # #bounds for bbp555
              # 0.0001 - pars[3],
              # pars[3] - 0.05,
              # 
              # #bounds for zB
              # 0.1 - pars[4],
              # pars[4] - 10,
              
              #bounds for aerial fractions
              pars[5:(length(pars))]
              #,
              
              # #bounds for zB
              # 0.0001 - pars[10],
              # pars[10] - 1
              
            )
            
          } else {
            fheq <- function(pars,data) sum(pars[5:(length(pars)-1)]) - 1
            
            fhin <- function(pars,data) c(#bounds for [chl]
              # 0.5 - pars[1],
              # pars[1] - 10,
              # 
              # #bounds for adg443
              # 0.1 - pars[2],
              # pars[2] - 10,
              # 
              # #bounds for bbp555
              # 0.0001 - pars[3],
              # pars[3] - 0.05,
              # 
              # #bounds for zB
              # 0.1 - pars[4],
              # pars[4] - 10,
              
              #bounds for aerial fractions
              pars[5:(length(pars)-1)]
              #,
              
              # #bounds for zB
              # 0.0001 - pars[10],
              # pars[10] - 1
              
            )
          }
          
          
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
          
          
          #print(MLE_estimates$par, digits=3)
          
          #print(sum(MLE_estimates$par[5:9]))
          
          #sprintf("%.5f", MLE_estimates$par)
          
        } else {
            
            MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata, 
                                  #lower = c(0, 0, 0),     # Lower bound on parameters
                                  #upper = c(10, 2, 0.3),  # Upper bound on parameters
                                  method = method.opt,
                                  control = list(parscale = abs(par0)),
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
      
      #@@@@@@@@@@@@@@ SSR Objective function @@@@@@@@@@@@@@@@@@@@@@@@
      if (obj.fn == "SSR") {
        print("Sum square of residual will be used to construct objective function")
        cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
        Fit <- data.frame("C_ph"=initial[1],
                          "a_cdom.440"=initial[2],
                          "a.nap.440"=initial[3])
        
        params <- as.numeric(Fit)
        cat(paste0("\033[0;36m","Initial values are ",params[1],"    ", params[2],"    ", 
                   params[3],"\033[0m","\n"))
        cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
        Sys.sleep(1)
        start.time = Sys.time()
        
        ##Single RUN
        inverse_output <- pracma::fminunc(fn = Saber_inverse,x0 = params,valdata=obsdata,
                                          bbp.550= Fit.input$bbp.550)
        
        cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
        
        Fit.optimized.ssobj <- data.frame("chl"=inverse_output$par[1], 
                                          "acdom.440"=inverse_output$par[2],
                                          "anap.440"=inverse_output$par[3])
        prmatrix(Fit.optimized.ssobj)
        
        #Calculate hessian matrix for var-covar matrix
        hessian.inverse <- numDeriv::hessian(x =as.numeric(Fit.optimized.ssobj), 
                                             func = Saber_inverse, valdata=obsdata,
                                             bbp.550= Fit.input$bbp.550)
        rownames(hessian.inverse) <- names(par0)
        colnames(hessian.inverse) <- names(par0)
        cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
        prmatrix(hessian.inverse)
        cat(paste0("\033[0;36m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
        
        MLE_par <- as.numeric(Fit.optimized.ssobj)
        MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
        MLE <- data.table("param" = names(par0),
                          "estimates" = MLE_par,
                          "sd(+/-)" = MLE_SE)
        
        print("The retrieved parameters are:")
        prmatrix(MLE)
        cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
        end.time = Sys.time(); time_taken <- end.time - start.time
        return(list(Fit.optimized.ssobj,"time.elapsed"= time_taken))
      }
    }
  }
  
  
  
}

#===============================================================================================================
#Solver for deep waters
#===============================================================================================================
solve.objective.inverse.deep.final <- function(
                                      bbp.constrain = F,
                                      bbp.550 =0.05,
                                               
                                     auto_spectral_slope = T,
                                     manual_spectral_slope = F,
                                     manual_spectral_slope_vals = c("s_g"=0.015, "s_d"=0.01160, "gamma"=0.5),

                                     obj.fn, initial, sa.model="am03", method.opt,
                                               
                                    obsdata,lower.b=c(0.1,0.1,0.0001), upper.b=c(30,10,1), 
                                    batch, pop.sd){
  if (auto_spectral_slope == TRUE & manual_spectral_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  
  myFun <- function(x) {
    NA
  }
  #initial <- abs(initial)
  if (obj.fn == "log-LL") {
    print("Log-Likelihood will be used to construct objective function")
    
    if (bbp.constrain == TRUE) {
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_final(
            use_true_IOPs = F, 
            a_non_water_path = IOP_files[idx_a],
            bb_non_water_path = IOP_files[idx_bb],
            
            chl = pars[1], 
            acdom440 =NULL, 
            anap440 =NULL , 
            a_dg = pars[2],
            bbp.550 = bbp.550,
            
            z = zB,
            use_spectral_rb = F, 
            rb.fraction = fA.set,
            
            
            realdata = data,
            
            slope.parametric = auto_spectral_slope,
            dg_composite = T,
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            use_spectral_shape_chl = F,
            use_spectral_shape_dg = T,
            
            sicf = F, q_phi = 0.02, 
            fDOM = F,
            verbose = F, plot = F
          )
          
        } else {
          Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                              anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
                              verbose = F, realdata = data, plot = F)
        }
        
        # Negative log-likelihood
        if(pop.sd == "TRUE") {
          if (batch == FALSE) {
            
            #ON for Albert Rrs
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) 
            
          } else {
            
            #ON for IOCCG
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE)) 
            
          }
        } else {
          if (batch == TRUE | batch == FALSE) {
            #browser()
            #ON for unknown sigma
            smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[3], log = TRUE)) 
            
          }
          
        }
        return(smull)
      }
    } else {
      ##Create log-likelihood function
      NLL_deep = function(pars, data) {
        
        # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward_final(
            use_true_IOPs = F, 
            a_non_water_path = IOP_files[idx_a],
            bb_non_water_path = IOP_files[idx_bb],
            
            chl = pars[1], 
            acdom440 =NULL, 
            anap440 =NULL , 
            a_dg = pars[2],
            bbp.550 = pars[3],
            
            z = zB,
            use_spectral_rb = F, 
            rb.fraction = fA.set,
            
            
            realdata = data,
            
            slope.parametric = auto_spectral_slope,
            dg_composite = T,
            
            use_manual_slope =manual_spectral_slope,
            manual_slope =  manual_spectral_slope_vals,
            
            use_spectral_shape_chl = F,
            use_spectral_shape_dg = T,
            
            sicf = F, q_phi = 0.02, 
            fDOM = F,
            verbose = F, plot = F
          )
          
        } else {
          Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                              anap440 =pars[3], bbp.550 = pars[3],
                              verbose = F, realdata = data, plot = F)
        }
        
        # Negative log-likelihood
        if(pop.sd == "TRUE") {
          if (batch == FALSE) {
            
            #ON for Albert Rrs
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) 
            
          } else {
            
            #ON for IOCCG
            smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE)) 
            
          }
        } else {
          if (batch == TRUE | batch == FALSE) {
            #browser()
            #ON for unknown sigma
            smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], log = TRUE)) 
            
          }
          
        }
        return(smull)
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
    
    #Turn ON known population sigma for likelihood
    if(pop.sd == TRUE) {
      if (bbp.constrain == TRUE) {
        par0 = c(chl = initial[1], adg440 = initial[2])#, bbp550 = initial[4]) 
      } else {
        par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3])
      }
      
    } else {
      #Turn ON for uknown population sigma for likelihood
      if (bbp.constrain == TRUE) {
        par0 = c(chl = initial[1], adg440 = initial[2], pop.sd = initial[3]) 
      } else {
        par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], pop.sd = initial[4])
      }
      
    }
    if (bbp.constrain == TRUE) {
      cat(paste0("\033[0;36m","Initial values are ",par0[1],"    ", par0[2],"\033[0m","\n"))
    } else {
      cat(paste0("\033[0;36m","Initial values are ",par0[1],"    ", par0[2],"    ", par0[3],"\033[0m","\n"))
    }
    
    
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
    
    if (obj.fn == "obj_L98") {
      print("Spectral error index from Lee et al. 1999 will be used to construct objective function")
      
      if (bbp.constrain == TRUE) {
        ##Create log-likelihood function
        NLL_deep = function(pars, data) {
          
          # Values predicted by the forward model for single RUN
          if (sa.model == "am03") {
            Gpred = Saber_forward_final(
              use_true_IOPs = F, 
              a_non_water_path = IOP_files[idx_a],
              bb_non_water_path = IOP_files[idx_bb],
              
              chl = pars[1], 
              acdom440 =NULL, 
              anap440 =NULL , 
              a_dg = pars[2],
              bbp.550 = bbp.550,
              
              z = zB,
              use_spectral_rb = F, 
              rb.fraction = fA.set,
              
              
              realdata = data,
              
              slope.parametric = auto_spectral_slope,
              dg_composite = T,
              
              use_manual_slope =manual_spectral_slope,
              manual_slope =  manual_spectral_slope_vals,
              
              use_spectral_shape_chl = F,
              use_spectral_shape_dg = T,
              
              sicf = F, q_phi = 0.02, 
              fDOM = F,
              verbose = F, plot = F
            )
            
          } else {
            Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                                anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
                                verbose = F, realdata = data, plot = F)
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
      } else {
        ##Create log-likelihood function
        NLL_deep = function(pars, data) {
          
          # Values predicted by the forward model for single RUN
          if (sa.model == "am03") {
            Gpred = Saber_forward_final(
              use_true_IOPs = F, 
              a_non_water_path = IOP_files[idx_a],
              bb_non_water_path = IOP_files[idx_bb],
              
              chl = pars[1], 
              acdom440 =NULL, 
              anap440 =NULL , 
              a_dg = pars[2],
              bbp.550 = pars[3],
              
              z = zB,
              use_spectral_rb = F, 
              rb.fraction = fA.set,
              
              
              realdata = data,
              
              slope.parametric = auto_spectral_slope,
              dg_composite = T,
              
              use_manual_slope =manual_spectral_slope,
              manual_slope =  manual_spectral_slope_vals,
              
              use_spectral_shape_chl = F,
              use_spectral_shape_dg = T,
              
              sicf = F, q_phi = 0.02, 
              fDOM = F,
              verbose = F, plot = F
            )
            
          } else {
            Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                                anap440 =pars[3], bbp.550 = pars[3],
                                verbose = F, realdata = data, plot = F)
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
      
      #Turn ON known population sigma for likelihood
      if(pop.sd == TRUE) {
        if (bbp.constrain == TRUE) {
          par0 = c(chl = initial[1], adg440 = initial[2])#, bbp550 = initial[4]) 
        } else {
          par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3])
        }
        
      } else {
        #Turn ON for uknown population sigma for likelihood
        if (bbp.constrain == TRUE) {
          par0 = c(chl = initial[1], adg440 = initial[2], pop.sd = initial[3]) 
        } else {
          par0 = c(chl = initial[1], adg440 = initial[2], bbp550 = initial[3], pop.sd = initial[4])
        }
        
      }
      if (bbp.constrain == TRUE) {
        cat(paste0("\033[0;36m","Initial values are ",par0[1],"    ", par0[2],"\033[0m","\n"))
      } else {
        cat(paste0("\033[0;36m","Initial values are ",par0[1],"    ", par0[2],"    ", par0[3],"\033[0m","\n"))
      }
      
      
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      #Sys.sleep(2)
      start.time = Sys.time()
      
      if (sa.model == "am03") {
        print("Albert & Mobley 2003 SA model used for Forward calculation")
      } else {
        print("Lee et al. 1999 SA model used for Forward calculation")
      }
      
      if (method.opt == "L-BFGS-B") {
        
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
    
    #@@@@@@@@@@@@@@ SSR Objective function @@@@@@@@@@@@@@@@@@@@@@@@
    if (obj.fn == "SSR") {
      print("Sum-square-of-residuals will be used to construct objective function")
      
      cat(paste0("\033[0;34m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
      Fit <- data.frame("C_ph"=initial[1],
                        "a_cdom.440"=initial[2],
                        "a.nap.440"=initial[3])
      
      params <- as.numeric(Fit)
      cat(paste0("\033[0;36m","Initial values are ",params[1],"    ", params[2],"    ", 
                 params[3],"\033[0m","\n"))
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      Sys.sleep(1)
      start.time = Sys.time()
      
      ##Single RUN
      inverse_output <- pracma::fminunc(fn = Saber_inverse,x0 = params,valdata=obsdata,
                                        bbp.550= Fit.input$bbp.550)
      cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
      Fit.optimized.ssobj <- data.frame("chl"=inverse_output$par[1], "acdom.440"=inverse_output$par[2],
                                        "anap.440"=inverse_output$par[3])
      prmatrix(Fit.optimized.ssobj)
      #Calculate hessian matrix for var-covar matrix
      hessian.inverse <- numDeriv::hessian(x =as.numeric(Fit.optimized.ssobj), 
                                           func = Saber_inverse, valdata=obsdata,
                                           bbp.550= Fit.input$bbp.550)
      rownames(hessian.inverse) <- names(par0)
      colnames(hessian.inverse) <- names(par0)
      cat(paste0("\033[0;32m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
      prmatrix(hessian.inverse)
      cat(paste0("\033[0;36m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
      
      MLE_par <- as.numeric(Fit.optimized.ssobj)
      MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
      MLE <- data.table("param" = names(par0),
                        "estimates" = MLE_par,
                        "sd(+/-)" = MLE_SE)
      
      print("The retrieved parameters are:")
      prmatrix(MLE)
      cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
      end.time = Sys.time(); time_taken <- end.time - start.time
      return(list(Fit.optimized.ssobj,"time.elapsed"= time_taken))
    }
  }
  
  
}
