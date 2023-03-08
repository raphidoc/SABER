#===============================================================================================
## inelastic.wavelength.redistribution.R contains wavelength redistribution functions (WRF) for
## SICF and fDOM and also wavelength discretization function for fDOM and SICF WRFs.

## The WRFs are coded folllowing Ocean Optics Web Book given formulations.
#===============================================================================================

#SICF wavelength redistribution function
wrfChl <- function(wave1, wave2, quant_phi=0.02) {
  sigmac <- 10.6
  wavec0 <- 685.0
  PhiChl = quant_phi
  factor1 <- 1.0/(sigmac*sqrt(2.0*pi))
  factor2 <- 0.5/(sigmac*sigmac)
  
  gchl <- ifelse(wave1 < 370.0 || wave1 > 690.0, 0.0, 1.0)
  
  hchl <- factor1 * exp(-factor2 * (wave2 - wavec0)^2)
  
  wrfchl <- PhiChl * gchl * hchl * wave1 / wave2
  
  return(wrfchl)
}


#fDOM wavelength redistribution function
wrfCDOM <- function(wave1, wave2, verbose_fdom = FALSE) {
  wave <- c(310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550)
  A0FA7 <- c(5.18e-5, 6.34e-5, 8.00e-5, 9.89e-5, 9.39e-5, 10.48e-5, 12.59e-5, 13.48e-5, 
             13.61e-5, 9.27e-5, 4.00e-5, 1.00e-5, 0.0)
  
  A1 <- 0.470
  A2 <- 0.407
  B1 <- 8.077e-4
  B2 <- -4.18e-4
  C1 <- -0.0181
  C2 <- -1.76
  
  if (wave1 >= wave[1] & wave1 <= wave[length(wave)]) {
    
    A0 = approx(x=wave, y=A0FA7, xout = wave1, method = "linear")$y
    temp = (1.0/wave2 - A1/wave1 - B1)/(0.6*(A2/wave1 + B2))
    wrfcdom = A0*exp(-temp*temp)*wave1/wave2
    
  } else {
    if (verbose_fdom == TRUE) {
      print(paste0("The current ex wavelength = ", wave1, " is not affected by the fDOM wavelength redistribution function"))
    }
    
    A0 = 0.0
    wrfcdom = 0.0
    
  }
  
  return(wrfcdom)
}

#Wavelength Discretization to perform the integration for the wavelength redistribution functions
wave.redist.inelastic <- function(wavelength = seq(400,800,10),
                                  include.chl = T,
                                  include.CDOM = T, quant_y_phi = 0.02) {
  
  deltaw <- 1
  waveb <- wavelength
  Nwave <- length(waveb) - 1
  
  if (include.chl == TRUE) {
    
    iChlfl <- 1
    print("SICF WRF is enabled for wavelength discretization")
  } else {
    
    iChlfl <- 0
  }
  
  if (include.CDOM == TRUE) {
    
    iCDOMfl <- 1
    print("fDOM WRF is enabled for wavelength discretization")
  } else {
    
    iCDOMfl <- 0
  }
  
  WRF_Chl <-matrix(0, nrow = Nwave, ncol = Nwave)
  WRF_CDOM <- matrix(0, nrow = Nwave, ncol = Nwave)
  WRF_Raman <- matrix(0, nrow = Nwave, ncol = Nwave)
  
  for (j in 2:Nwave) {
    wj <- waveb[j] + 0.5 * deltaw
    nj <- as.integer((waveb[j + 1] - waveb[j])/deltaw)
    factor <- deltaw * deltaw / (waveb[j + 1] - waveb[j])
    
    for (i in 1:(j - 1)) {
      wi <- waveb[i] + 0.5 * deltaw
      ni <- as.integer((waveb[i + 1] - waveb[i])/deltaw)
      
      if (iChlfl != 0) {
        sumchl <- 0
        for (jj in 1:nj) {
          wavej <- wj + as.numeric(jj - 1) * deltaw
          for (ii in 1:ni) {
            wavei <- wi + as.numeric(ii - 1) * deltaw
            sumchl <- sumchl + wrfChl(wavei, wavej, quant_phi = quant_y_phi)
          }
        }
        
        WRF_Chl[i, j] <- factor * sumchl
      }
      
      if (iCDOMfl != 0) {
        sumcdom <- 0
        for (jj in 1:nj) {
          wavej <- wj + as.numeric(jj - 1) * deltaw
          for (ii in 1:ni) {
            wavei <- wi + as.numeric(ii - 1) * deltaw
            sumcdom <- sumcdom + wrfCDOM(wavei, wavej)
          }
        }
        
        WRF_CDOM[i, j] <- factor * sumcdom
      }
    }
  }
  return(list("sicf"=WRF_Chl, "fDOM"=WRF_CDOM))
}




