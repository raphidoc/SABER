#' inversion_unconstrained
#' Unconstrained full inversion. Two options for forward semi analytical model
#'  and two option for the objective function.
#'
#' @author Soham Mukherjee
#'
#' @param sa_model c("am03", "lee98")
#' @param objective_fct c("log-ll", "SSR", "lee99") SSR is not implemented in the code
#' @param method_opt c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent","levenberg-marqardt", "auglag")
#'
#' @export


inversion_unconstrained <- function(
    wavelength,
    sa_model,
    objective_fct
    ) {

  rlang::inform(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  rlang::inform(paste0("\033[0;39m","########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######","\033[0m","\n"))
  rlang::inform(paste0("\033[0;32m","###################################################################","\033[0m","\n"))

# Estimate Rrs from a set of starting parameter ---------------------------

  if (sa_model == "am03") {
    init_rrs <- saber_forward(
      use_true_IOPs = F,
      chl = pars[1],
      a_dg = pars[2],
      bbp.550 = pars[3],
      z = pars[4],
      rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
      Rrs_input_for_slope = data,
      slope.parametric = auto_spectral_slope,
      use_manual_slope =manual_spectral_slope,
      manual_slope =  manual_spectral_slope_vals,
      verbose = F, wavelength = wavelength
    )
  } else if (sa_model == "lee98") {
    init_rrs <- lee98_forward(
      wavelength = wavelength,

    )

  } else {
    rlang::abort("sa_model must be either am03 or lee98")
  }

  if (objective_fct == "log-ll") {

    # Negative log-likelihood
    smull <- -sum(
      dnorm(
        x = 10000*data,
        mean = 10000*init_rrs$rrs_0m,
        sd = pars[length(pars)], ###### ????????
        log = TRUE)
      )

  } else if (objective_fct == "lee99") {

    #The Spectral error index from Lee 1999
    # Define the spectral regions
    region1 <- which(wavelength >= 400 & wavelength <= 675)
    region2 <- which(wavelength >= 750 & wavelength <= 830)

    # Calculate the numerator of the error index
    numerator <- sqrt(sum((data[region1] - init_rrs$rrs_0m[region1])^2) +
                        sum((data[region2] - init_rrs$rrs_0m[region2])^2))

    # Calculate the denominator of the error index
    denominator <- sum(init_rrs$rrs_0m[region1]) + sum(init_rrs$rrs_0m[region2])

    # Calculate the error index
    smull <- numerator / denominator
  } else {
    rlang::abort("objective function must be either 'log-ll' or 'lee99' ")
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
  par0 = c(
    chl = initial[1],
    adg440 = initial[2],
    bbp550 = initial[3],
    z = initial[4],
    initial[5:(4+initial_rb_length)],
    pop.sd = initial[length(initial)]
  )


  # Optimization ------------------------------------------------------------


  start.time = Sys.time()

  if (sa_model == "am03") {
    print("Albert & Mobley 2003 SA model used for forward model")
  } else {
    print("Lee et al. 1999 SA model used for forward model")
  }

  if (method_opt == "L-BFGS-B") {

    MLE_estimates = optim(par = par0, fn = NLL_unconstr, data = obsdata,
                          lower = lower_b,     # Lower bound on parameters
                          upper = upper_b,  # Upper bound on parameters
                          #gr = grad.calculate,
                          method = method_opt,
                          control = list(#fnscale = 1,
                            parscale = abs(par0)),
                          #hessian = TRUE
    )
  }
  if (method_opt == "levenberg-marqardt") {

    LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL_unconstr,data = obsdata, print.info = F)
    MLE_estimates = data.frame("par"=LM_estimates$b)

  }

  if (method_opt == "auglag") {
    print("Augmented Lagriangian with equality constraints will be used for inversion")


    fheq <- function(pars,data) sum(pars[5:(length(pars)-1)]) - 1

    fhin <- function(pars,data) c(
      #bounds for areal fractions
      pars[5:(length(pars)-1)]
    )

    MLE_estimates <- alabama::auglag(
      fn = NLL_unconstr,
      par = par0,
      heq = fheq,
      hin = fhin,
      data = obsdata,
      control.outer = list(trace = F, method = "nlminb")
    )
  }

  if (method_opt == "Nelder-Mead" | method_opt ==  "SANN" | method_opt ==  "Brent") {

    print("Simplex based Optimization will be used for inversion")
    MLE_estimates = optim(
      par = par0,
      fn = NLL_unconstr,
      data = obsdata,
      method = method_opt,
      control = list(
        parscale = as.numeric(abs(par0))
      ),
      hessian = FALSE
    )
  }

# Calculate uncertainty ---------------------------------------------------

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

  end.time = Sys.time()

  if (!is.numeric(result)) {
    cat(paste0("\033[0;31m","Failed to calculate diagonal of hessian from high degree of correlation, coerce to NA","\033[0m","\n"))
    MLE_SE = result
  } else {
    MLE_SE <- result
  }
  MLE <- data.table(
    "param" = names(par0),
    "estimates" = MLE_par,
    "sd(+/-)" = MLE_SE)

  if (verbose) {
    rlang::inform("The retrieved parameters are:")
    prmatrix(MLE)

    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
      rlang::inform(paste0("\033[0;32m","CONVERGENCE: GLOBAL","\033[0m","\n"))
    } else {
      convergence = "FALSE"
      rlang::inform(paste0("\033[0;34m","CONVERGENCE: LOCAL","\033[0m","\n"))
    }

    rlang::inform(glue:glue(("time.elapsed: ", time_taken))


    time_taken <- end.time - start.time
    return(list(MLE, "convergence"= convergence))
  }

  cat(paste0("\033[0;32m","#################### INVERSION ENDS #########################","\033[0m","\n"))
}
