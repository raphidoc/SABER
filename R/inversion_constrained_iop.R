# IOP constrained inversion -----------------------------------------------

if (constrain_shallow_iop == TRUE) {
  initial_rb_length = length(initial[2:(length(initial)-1)])

  if (obj_fn == "log-LL") {
    print("Log-Likelihood will be used to construct objective function")
    ##Create log-likelihood function

    NLL_unconstr = function(pars, data) {

      # Values predicted by the forward model for single RUN
      if (sa_model == "am03") {
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

          verbose = F, wavelength = wavelength,plot = F
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

          verbose = F, wavelength = wavelength, plot = F,
          type_Rrs_below = "shallow"
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
      if (sa_model == "am03") {
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

          verbose = F, wavelength = wavelength, plot = F
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
          verbose = F, wavelength = wavelength, plot = F,
          type_Rrs_below = "shallow"
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

  if (sa_model == "am03") {
    print("Albert & Mobley 2003 SA model used for forward model")
  } else {
    print("Lee et al. 1999 SA model used for forward model")
  }

  if (method_opt == "L-BFGS-B") {

    print("L-BFGS-B Optimization will be used for inversion")
    MLE_estimates = optim(par = as.numeric(par0), fn = NLL_unconstr, data = obsdata,
                          lower = as.numeric(lower_b),     # Lower bound on parameters
                          upper = as.numeric(upper_b),  # Upper bound on parameters
                          #gr = grad.calculate,
                          method = method_opt,
                          control = list(#fnscale = 1,
                            parscale = as.numeric(abs(par0))),
                          #hessian = TRUE
    )
  }

  if (method_opt == "levenberg-marqardt") {

    print("Levenberg-Marquardt Optimization will be used for inversion")
    LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_unconstr,data = obsdata, print.info = F)
    MLE_estimates = data.frame("par"=LM_estimates$b)
  }

  if (method_opt == "auglag") {
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
                                     #lower = lower_bound,
                                     #upper = upper_bound,
                                     #data=surface_rrs_translate(insitu.data),
                                     data = obsdata,
                                     control.outer = list(trace = F, method = "nlminb")
    )





    #print(inverse_output_fmincon$par, digits=3)

    #print(sum(inverse_output_fmincon$par[5:9]))

    #sprintf("%.5f", inverse_output_fmincon$par)

  }

  if (method_opt == "Nelder-Mead" |method_opt ==  "SANN" | method_opt ==  "Brent") {

    print("Simplex based Optimization will be used for inversion")
    MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata,
                          #lower = c(0, 0, 0),     # Lower bound on parameters
                          #upper = c(10, 2, 0.3),  # Upper bound on parameters
                          method = method_opt,
                          control = list(#fnscale = 1,
                            parscale = as.numeric(abs(par0_constrained))),
                          hessian = FALSE)


  }


  cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))

  #print(MLE_estimates$par)

  cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
  Sys.sleep(1)

  #Calculate hessian matrix for var-covar matrix
  if (method_opt == "auglag") {
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
