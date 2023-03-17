#=============================================================================
# grad.code.R consists of code to calculate the gradient function analytically
#for the likelihood function optimization 
#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=============================================================================

NLL = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  
      Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                            anap440 =pars[3], bbp.550 =Fit.input$bbp.550,
                            verbose = F, realdata = data, plot = F)
      
  # # # Values predicted by the forward model for BATCH RUN
  #     Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
  #                           anap440 =pars[3], bbp.550 = HL.deep.iop$bbp550[j],
  #                           verbose = F, realdata = data, plot = F)
  
  
  # Negative log-likelihood 

   ll= -sum(dnorm(x = 1000*data, mean = 1000*Gpred[[1]]$Rrs, sd = pars[4], log = TRUE)) #ON for unknown sigma
   #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
   #attr(ll, "gradient") <- grad
   return(ll)
  
}

#par0 = par0[-3]
lower.bound <- par0 - 0.8*par0
upper.bound <- par0 + 5*par0

L_BFS_estimates = optim(par = par0, fn = NLL, data = obsdata, 
                      lower = lower.bound,     # Lower bound on parameters
                      upper = upper.bound,  # Upper bound on parameters
                      gr = grad_NLL,
                      
                      method = methods.opt[4],
                      control = list(parscale = abs(par0)),
                      hessian = F)

LM_estimates = marqLevAlg::marqLevAlg(b = par0, fn = NLL,data = obsdata, print.info = F)

LM_estimates = marqLevAlg::marqLevAlg(b = c(5,1,0.001), fn = NLL,data = obsdata, print.info = T)
summary(LM_estimates, log=T)

numDeriv::grad(x=par0*0.2, func=NLL,data= obsdata)

#Try optimization with nloptr
NLL_nloptr_mle = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  
  Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
                        z = pars[4], 
                        # rb0 = pars[5],
                        # rb1 = pars[6],
                        # rb2 = pars[7],
                        # rb3 = pars[8],
                        # rb4 = pars[9],
                        # rb5 = pars[10],
                        
                        rb.fraction = as.numeric(pars[5:(length(pars)-2)]),
                        verbose = F, realdata = data, plot = F)
  
  # # # Values predicted by the forward model for BATCH RUN
  #     Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
  #                           anap440 =pars[3], bbp.550 = HL.deep.iop$bbp550[j],
  #                           verbose = F, realdata = data, plot = F)
  
  
  # Negative log-likelihood 
  
  ll= -sum(dnorm(x = 1000*data, mean = 1000*Gpred[[1]]$Rrs, sd = pars[length(pars)], log = TRUE)) #ON for unknown sigma
  #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
  #attr(ll, "gradient") <- grad
  return(ll)
  
}

#Try optimization with nloptr
NLL_nloptr_sse = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  
  # Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
  #                       anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
  #                       z = pars[4], 
  #                       # rb0 = pars[5],
  #                       # rb1 = pars[6],
  #                       # rb2 = pars[7],
  #                       # rb3 = pars[8],
  #                       # rb4 = pars[9],
  #                       # rb5 = pars[10],
  #                       
  #                       rb.fraction = as.numeric(pars[5:(length(pars)-2)]),
  #                       verbose = F, realdata = data, plot = F)
  
  Gpred = invisible(Saber_forward_paramteric_conc(chl = pars[1], 
                                        acdom440 =NULL, 
                                        anap440 =NULL , 
                                        a_dg = pars[2],
                                        bbp.550 = Fit.input$bbp.550,
                                        z = pars[3],
                                        #realdata = rrs.forward.am,
                                        #realdata = insitu.data,
                                        dg_composite = T,
                                        slope.parametric = T,
                                        use_spectral_shape_chl = F,
                                        use_spectral_shape_dg = T,
                                        sicf = F, q_phi = 0.02,
                                        rb.fraction = as.numeric(pars[4:(length(pars))]),
                                        verbose = F, realdata = data, plot = F))
  

  
  rrs_est = Gpred[[1]]$Rrs
  # Sum-squared residual of error (SSE)
  
  sse= sum((data - rrs_est)^2)
  #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
  #attr(ll, "gradient") <- grad
  return(sse)
  
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

#Constrained
#Try optimization with nloptr
NLL_nloptr_sse_constr = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  
  # Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
  #                       anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
  #                       z = pars[4], 
  #                       # rb0 = pars[5],
  #                       # rb1 = pars[6],
  #                       # rb2 = pars[7],
  #                       # rb3 = pars[8],
  #                       # rb4 = pars[9],
  #                       # rb5 = pars[10],
  #                       
  #                       rb.fraction = as.numeric(pars[5:(length(pars)-2)]),
  #                       verbose = F, realdata = data, plot = F)
  
  Gpred = invisible(Saber_forward_paramteric_conc(chl = Fit.input$chl, 
                                        acdom440 =Fit.input$acdom.440, 
                                        anap440 = Fit.input$anap.440 , 
                                        a_dg = NULL,
                                        bbp.550 = Fit.input$bbp.550,
                                        z = pars[1],
                                        #realdata = rrs.forward.am,
                                        #realdata = insitu.data,
                                        dg_composite = F,
                                        slope.parametric = F,
                                        use_spectral_shape_chl = F,
                                        use_spectral_shape_dg = F,
                                        sicf = F, q_phi = 0.02,
                                        rb.fraction = as.numeric(pars[2:(length(pars))]),
                                        verbose = F, realdata = data, plot = F))
  
  
  
  rrs_est = Gpred[[1]]$Rrs
  # Sum-squared residual of error (SSE)
  
  sse= sum((data - rrs_est)^2)
  #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
  #attr(ll, "gradient") <- grad
  return(sse)
  
}

numDeriv::grad(x=par0_constr*0.2, func=NLL_nloptr_sse_constr,data= obsdata)

grad_sse <- function(pars, data) {
  
  # Predicted Rrs values
  Gpred <- Saber_forward_paramteric_conc(chl = Fit.input$chl, 
                                         acdom440 =Fit.input$acdom.440, 
                                         anap440 = Fit.input$anap.440 , 
                                         a_dg = NULL,
                                         bbp.550 = Fit.input$bbp.550,
                                         z = pars[1],
                                         rb.fraction = as.numeric(pars[2:(length(pars))]),
                                         verbose = F, realdata = data, plot = F)
  rrs_est <- Gpred[[1]]$Rrs
  
  # Gradient with respect to model parameter "z"
  dz <- -2 * sum(data - rrs_est) * Gpred[[2]][1]
  
  # Gradient with respect to model parameter "rb.fraction"
  drbf <- -2 * t(Gpred[[1]]$dR_drb) %*% (data - rrs_est)
  drbf <- c(drbf[1], drbf[, 2:length(drbf)] %*% diag(as.numeric(pars[2:(length(pars))] > 0)))
  
  # Combine gradients into single vector
  grad <- c(dz, drbf)
  return(grad)
}


eval_grad_f <- function(pars,data) {
  Gpred = Saber_forward_paramteric_conc(chl = pars[1], 
                                        acdom440 =NULL, 
                                        anap440 =NULL , 
                                        a_dg = pars[2],
                                        bbp.550 = Fit.input$bbp.550,
                                        z = pars[3],
                                        #realdata = rrs.forward.am,
                                        #realdata = insitu.data,
                                        dg_composite = T,
                                        slope.parametric = T,
                                        use_spectral_shape_chl = F,
                                        use_spectral_shape_dg = T,
                                        sicf = F, q_phi = 0.02,
                                        rb.fraction = as.numeric(pars[4:(length(pars))]),
                                        verbose = F, realdata = data, plot = F)
  
  
  
  rrs_est = Gpred[[1]]$Rrs
  
  # Sum-squared residual of error (SSE)
  #sse_gr = Deriv(substitute(sum((data - rrs_est)^2)),x = c("data", "rrs_est"))
  grad_rrs_est <- 2 * (rrs_est - data)
  grad_data <- -2 * (rrs_est - data)
  #return(c(grad_rrs_est, grad_data))
  return(c(grad_rrs_est))
}

eval_grad_f_constr <- function(pars,data) {
  Gpred = invisible(Saber_forward_paramteric_conc(chl = Fit.input$chl, 
                                                  acdom440 =Fit.input$acdom.440, 
                                                  anap440 = Fit.input$anap.440 , 
                                                  a_dg = NULL,
                                                  bbp.550 = Fit.input$bbp.550,
                                                  z = pars[1],
                                                  #realdata = rrs.forward.am,
                                                  #realdata = insitu.data,
                                                  dg_composite = F,
                                                  slope.parametric = F,
                                                  use_spectral_shape_chl = F,
                                                  use_spectral_shape_dg = F,
                                                  sicf = F, q_phi = 0.02,
                                                  rb.fraction = as.numeric(pars[2:(length(pars))]),
                                                  verbose = F, realdata = data, plot = F))
  
  
  
  rrs_est = Gpred[[1]]$Rrs
  
  # Sum-squared residual of error (SSE)
  #sse_gr = Deriv(substitute(sum((data - rrs_est)^2)),x = c("data", "rrs_est"))
  grad_rrs_est <- 2 * (rrs_est - data)
  grad_data <- -2 * (rrs_est - data)
  #return(c(grad_rrs_est, grad_data))
  return(c(grad_rrs_est))
}


# Equality constraints
eval_g_eq <- function(pars,data)
{
  return ( pars[4] + pars[5] +pars[6] + pars[7] + pars[8] - 1 )
}

# Equality constraints on constrained inversion
eval_g_eq_constr <- function(pars,data)
{
  return ( pars[2] + pars[3] +pars[4] + pars[5] + pars[6] - 1 )
}

# Lower and upper bounds
par0 = c(chl = 8, 
         a_dg = 1.5,
         #acdom440 = 0.8,
         #anap440 = 0.05,
         z = 2.5, 
         #rb.0 = 0.1,
         rb.1 = 0.5,
         rb.2 = 0.5,
         rb.3 = 0.5,
         rb.4 = 0.5,
         rb.5 = 0.5)
         #,population.sd = 0.1)

par0_constr = c(
         z = 2.5, 
         #rb.0 = 0.1,
         rb.1 = 0.5,
         rb.2 = 0.5,
         rb.3 = 0.5,
         rb.4 = 0.5,
         rb.5 = 0.5)
#,population.sd = 0.1)

#Autoscale Intital values from pre-Ft
increament.scale <- 1

lower.bound <- c((par0[1:2] - 0.8*par0[1:2]),z = 0.1,
                 #rb.0 = 0.1,
                 rb.1 = 0,
                 rb.2 = 0,
                 rb.3 = 0,
                 rb.4 = 0,
                 rb.5 = 0)
                # ,population.sd = 0.0001)

upper.bound <- c((par0[1:2] + 5*par0[1:2]),z = 10,
                 #rb.0 = 1,
                 rb.1 = 1,
                 rb.2 = 1,
                 rb.3 = 1,
                 rb.4 = 1,
                 rb.5 = 1)
                 #,population.sd = 1)

lower.bound_constr <- c(z = 1,
                 #rb.0 = 0.1,
                 rb.1 = 0,
                 rb.2 = 0,
                 rb.3 = 0,
                 rb.4 = 0,
                 rb.5 = 0)
# ,population.sd = 0.0001)

upper.bound_constr <- c(z = 10,
                 #rb.0 = 1,
                 rb.1 = 1,
                 rb.2 = 1,
                 rb.3 = 1,
                 rb.4 = 1,
                 rb.5 = 1)
#,population.sd = 1)


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-5 )

opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-5,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )

#optimization with gradient
opts <- list("algorithm"="NLOPT_LD_LBFGS",
             "xtol_rel"=1.0e-8)

res <- nloptr::nloptr ( x0 = par0,
                eval_f = NLL_nloptr_sse,
                #eval_grad_f = eval_grad_f,
                data = obsdata,
                lb = lower.bound,
                ub = upper.bound,
                #eval_g_ineq = eval_g_ineq,
                eval_g_eq = eval_g_eq,
                opts = opts
)
print(res)

res_constr <- nloptr::nloptr ( x0 = par0_constr,
                        eval_f = NLL_nloptr_sse_constr,
                        eval_grad_f = eval_grad_f_constr,
                        data = obsdata,
                        lb = lower.bound_constr,
                        ub = upper.bound_constr,
                        #eval_g_ineq = eval_g_ineq_constr,
                        eval_g_eq = eval_g_eq_constr,
                        opts = opts
)
print(res_constr)

solution_nloptr = data.frame(t(res$solution))
names(solution_nloptr) = c("chl", "a_dg", "z", "rb1", "rb2", "rb3", "rb4", "rb5")

saber_forward_param_nloptr = Saber_forward_paramteric_conc(chl = solution_nloptr$chl, 
                              acdom440 =NULL, 
                              anap440 =NULL , 
                              a_dg = solution_nloptr$a_dg,
                              bbp.550 = Fit.input$bbp.550,
                              z =solution_nloptr$z,
                              realdata = rrs.forward.am,
                              #realdata = insitu.data,
                              slope.parametric = T,
                              dg_composite = T,
                              use_spectral_shape_chl = F,
                              use_spectral_shape_dg = T,
                              sicf = F, q_phi = 0.02,
                              rb.fraction = as.numeric(solution_nloptr[4:(length(solution_nloptr))]),
                              verbose = F, plot = F)
plot(wavelength, obsdata)
lines(wavelength, saber_forward_param_nloptr[[1]]$Rrs)
#------------------------------------------------------------------------------------

saber.f.grad = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  #browser()
  Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 =Fit.input$bbp.550,
                        verbose = F, realdata = data, plot = F)
  
  # # # Values predicted by the forward model for BATCH RUN
  #     Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
  #                           anap440 =pars[3], bbp.550 = HL.deep.iop$bbp550[j],
  #                           verbose = F, realdata = data, plot = F)
  
  
  # Negative log-likelihood 
  
  ll= Gpred[[1]]$Rrs
  #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
  #attr(ll, "gradient") <- grad
  return(sum(ll))
  
}

grad.calculate <- function(data,par) {
  
  gradlist = numDeriv::grad(x=par, func=saber.f.grad, data= data)
  
  return(data.frame("dRrs_dchl" = gradlist[1],
                    "dRrs_dacdom440 " = gradlist[2],
                    "dRrs_danap440" = gradlist[3]))
}

test_grad = grad.calculate(par = c(4.77, 0.92, 0.01), data = obsdata)

numDeriv::grad(x=c(4.77,0.92,0.01), func=saber.f.grad, data= obsdata )

grad_NLL <- function(pars, data) {
  # Forward model predictions
  Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                      anap440 =pars[3], bbp.550 =Fit.input$bbp.550,
                      verbose = F, realdata = data, plot = F)
  
  test_grad = grad.calculate(par = pars[-length(pars)], data = obsdata)
  test_grad
  
  # Gradients for each parameter
  grad_chl = -sum((data - Gpred[[1]]$Rrs) * test_grad$dRrs_dchl) / pars[4]^2
  grad_acdom440 = -sum((data - Gpred[[1]]$Rrs) * test_grad$dRrs_dacdom440) / pars[4]^2
  grad_anap440 = -sum((data - Gpred[[1]]$Rrs) * test_grad$dRrs_danap440) / pars[4]^2
  grad_sigma = -sum((data - Gpred[[1]]$Rrs) * (data - Gpred[[1]]$Rrs) / pars[4]^3 - 1 / pars[4])
  
  # Return the gradient as a list
  return(data.frame("grad_chl"=grad_chl, 
                    "grad_cdom"=grad_acdom440, 
                    "grad_nap"=grad_anap440, 
                    "grad_sigma"=grad_sigma, row.names = "grad"))
}

ll_grad = grad_NLL(pars = par0, data = obsdata)


numDeriv::grad(x=par0, func=NLL,data= obsdata)



# diff(dnorm(x = obsdata, mean = rrs.forward.am, sd = 0.001, log = TRUE)) / diff(x = obsdata)
# 
# ddnorm <- function(x) eval(DD(expression(-sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, 
#                                                     sd = pars[4], log = TRUE))), "x", 
#                               order = 1))
# ddnorm(x=x)
# 
# mlogl4 <- function(theta, x) {
#    if (length(theta) != 2)
#      stop("length(theta) must be 2")
#    mu <- theta[1]
#    sigma <- theta[2]
#    value <- sum(-dcauchy(x, location = mu, scale = sigma, log = TRUE))
#    denom <- sigma^2 + (x - mu)^2
#    grad1 <- sum(-2 * (x - mu)/denom)
#    grad2 <- sum(-(1/sigma - (2 * sigma/denom)))
#    attr(value, "gradient") <- c(grad1, grad2)
#    return(value)
# }
# 
# out <- mlogl4(c(0.5,0.1),x = x)
# 
# mlogl4 <- function(theta, x) {
#   if (length(theta) != 2)
#     stop("length(theta) must be 2")
#   mu <- theta[1]
#   sigma <- theta[2]
#   value <- sum(-dcauchy(x, location = mu, scale = sigma, log = TRUE))
#   denom <- sigma^2 + (x - mu)^2
#   grad1 <- sum(-2 * (x - mu)/denom)
#   grad2 <- sum(-(1/sigma - (2 * sigma/denom)))
#   attr(value, "gradient") <- c(grad1, grad2)
#   return(value)
# }
# 
# out <- mlogl4(c(0.5,0.1),x = x)




# Function to calculate the gradient of the NLL function for a normal distribution
NLL_gradient <- function(params, obsdata, Saber_forward) {
  # Extract parameters
  mu <- Saber_forward(params)
  sigma <- params[length(params)]
  
  # Calculate NLL
  NLL <- sum(dnorm(obsdata, mean = mu, sd = sigma, log = TRUE))
  
  # Calculate gradient of NLL with respect to each parameter
  grad_mu <- sum((obsdata - mu) / (sigma^2))
  grad_sigma <- sum(((obsdata - mu)^2 - sigma^2) / (sigma^3))
  grad_params <- c(grad_mu, grad_sigma)
  
  # Return gradient of NLL
  return(grad_params)
}

Saber_gradient <- function(obsdata, params, sigma) {
  # calculate mean from the Saber_forward function
  mean <- Saber_forward(chl = params[1],
                        acdom440 = params[2],
                        anap440 = params[3],
                        bbp.550 = Fit.input$bbp.550)
  mean = mean[[1]]$Rrs
  # calculate the gradient of the log-likelihood function
  dNLL_dmean <- sum((obsdata - mean) / (sigma^2))
  dNLL_dsigma <- (length(obsdata) + 1) / sigma + sum((obsdata - mean)^2) / (sigma^3)
  
  # return the gradient
  return(c(dNLL_dmean, dNLL_dsigma))
}

Saber_gradient(obsdata = obsdata, params = c(10,5,1), sigma = 0.01)


# Function to calculate gradient of negative log-likelihood
nll_gradient <- function(pars, data) {
  mean <- Saber_forward(chl = pars[1],
                      acdom440 = pars[2],
                      anap440 = pars[3],
                      bbp.550 = Fit.input$bbp.550)
  mu = mean[[1]]$Rrs
  sigma <- params[length(params)]
  
  # Calculate the negative log-likelihood
  nll <- -sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
  
  # Calculate gradient with respect to mu
  grad_mu <- sum((data - mu) / (sigma^2)) * exp(-nll)
  
  # Calculate gradient with respect to sigma
  grad_sigma <- sum(((data - mu)^2 - sigma^2) / (sigma^3)) * exp(-nll)
  return(c(grad_mu, grad_sigma))
}

nll_gradient(pars = c(10,5,1,0.01), data = obs.data)

#====================================================================================


inverse_output <- pracma::fminunc(x0 = par0, fn = NLL,data=obsdata)


