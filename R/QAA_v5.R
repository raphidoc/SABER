QAA.v5 <- function(waves, Rrs, translate_subsurf=T)
{
  #wl  <- c(412,443,490,510,555,670)
  bbw <- Riops::spectral.bw(waves) * 0.5
  aw <- Riops::spectral.aw(waves)
  
  ix411 = which.min(abs(waves-411))
  ix443 = which.min(abs(waves-443))
  ix490 = which.min(abs(waves-490))
  ix510 = which.min(abs(waves-510))
  ix555 = which.min(abs(waves-555))
  ix640 = which.min(abs(waves-640))
  ix667 = which.min(abs(waves-667))
  
  nwaves = length(waves)
  a <- rep(0,nwaves)
  bb <- rep(0,nwaves)
  a2 <- rep(0,nwaves)
  bb2 <- rep(0,nwaves)
  g0 = 0.0895
  g1 = 0.1247
  
  # Check if 667 chanel is available
  if (min(abs(waves-667)) > 5 ) {
    print("No 667 nm channel")
    Rrs667 = 1.27*Rrs[ix555]^1.47 + 0.00018*(Rrs[ix490]/Rrs[ix555])^-3.19
    rrs667 = Rrs667 / (0.52 + 1.7*Rrs667)
  } else
  {
    up.lim <- 20 * Rrs[ix555]^1.5
    low.lim <- 0.9 * Rrs[ix555]^1.7
    if (Rrs[ix667] > up.lim | Rrs[ix667] < low.lim ) {
      Rrs667 = 1.27*Rrs[ix555]^1.47 + 0.00018*(Rrs[ix490]/Rrs[ix555])^-3.19
      rrs667 = Rrs667 / (0.52 + 1.7*Rrs667)
    } else {
      Rrs667 = Rrs[ix667]
      rrs667 = Rrs[ix667] / (0.52 + 1.7*Rrs[ix667])
    }
  }
  
  if (translate_subsurf == T) {
    #  STEP 0 - compute RRS to below Sea Surface
    rrs <-  Rrs / (0.52 + 1.7*Rrs)
    
  } else {
    rrs = Rrs
  }
  
  
  # Step 1 - Compute  bb/a+bb ratio
  X <-  (-g0 + sqrt(g0^2 + 4*g1*rrs)) / (2*g1)
  
  # Step 2 - Empirical estimation of a(555) from v5
  xsi <- log10((rrs[ix443]+rrs[ix490])/
                 (rrs[ix555]+5*rrs667/rrs[ix490]*rrs667))
  a.ref <- aw[ix555] + 10^(-1.146 - 1.366*xsi -0.469*xsi^2)
  
  # Step 3  - Analytical estimation of bb at 555 nm
  bbp.ref <- (X[ix555]*a.ref)/(1-X[ix555]) - bbw[ix555]
  
  # Step 4 - Empirical estimation of Nu
  nu <- 2.0*(1 - 1.2*exp(-0.9*rrs[ix443]/rrs[ix555]))
  
  # Step 5 - Spectral bbp
  bb <- bbw + bbp.ref*(waves[ix555]/waves)^nu
  
  # Step 6 - Spectral a
  a <-  ((1 - X)*bb)/X
  
  #parametric formula to retrieve component specific absorption
  s_cdm = 0.015 + (0.002/(0.6 + (Rrs[which.min(abs(wavelength - 443))]/Rrs[which.min(abs(wavelength - 555))])))
  epsilon_prime = exp(s_cdm*(443-411))
  
  tau_prime = 0.74 + (0.2/(0.8+(Rrs[which.min(abs(wavelength - 443))]/Rrs[which.min(abs(wavelength - 555))])))
  
  
  a_dg_443 = ((a[ix411]-(tau_prime*a[ix443])) - (aw[ix411]-(tau_prime*aw[ix443])))/(epsilon_prime - tau_prime)
  
  a_phi = a - aw - (a_dg_443*exp(-s_cdm*(wavelength - 443)))
  a_phi_443 = a_phi[ix443]
  
  chl_est = (a_phi_443/0.06)^(1/0.65) #Estimate Chlorophyll-a
  bbp_est = bb - bbw #Estimate spectral bbp
  a_nw_est = a - aw #Estimate spectral non-water absorption
  
  return(list(a=a,bb=bb, bbp = bbp_est, a_nw = a_nw_est,
              s_cdm=s_cdm, a_dg_443 = a_dg_443, 
              a_phi= a_phi, chl = chl_est,
              b_bp_555=bbp.ref))
  
}

# qaa_output = QAA.v5(waves = wavelength, Rrs = insitu.data)
# 
# #parametric formula to retrieve spectral slope of CDOM + NAP
# s_cdm = 0.015 + (0.002/(0.6 + (Rrs[which.min(abs(wavelength - 443))]/Rrs[which.min(abs(wavelength - 555))])))
# epsilon_prime = exp(s_cdm*(443-411))
# 
# tau_prime = 0.74 + (0.2/(0.8+(Rrs[which.min(abs(wavelength - 443))]/Rrs[which.min(abs(wavelength - 555))])))
# 
# 
# a_dg_443 = ((a[ix411]-(tau_prime*a[ix443])) - (aw[ix411]-(tau_prime*aw[ix443])))/(epsilon_prime - tau_prime)
# 
# a_phi = a - aw - (a_dg_443*exp(-s_cdm*(wavelength - 443)))
#   
#   
  
  
  
  