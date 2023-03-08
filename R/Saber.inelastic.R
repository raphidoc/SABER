
Rrs_Fluorescence <- function(dg_comsposite=TRUE, dg_443,
                             c_chl, abs_cdom_443, abs_nap_443,
                             Ed_path = "./input-spectra/Ed_HL.csv",
                             use_analytic_Ed = TRUE,
                             coeff_x=c(0.0992,0.40,0.078), 
                             wavelength=seq(400,800,10), phi_f=0.01) {
  
  ##approximated model :  A. Gilerson, J. Zhou, S. Hlaing, I. Ioannou, J. Schalles, 
  #B. Gross, F. Moshary, and S. Ahmed, "Fluorescence component in the reflectance spectra 
  #from coastal waters. Dependence on water composition," Opt. Express 15, 
  #15702-15721 (2007)
  if (Ed_analytical == TRUE) {
    
  }
  Ed_sim = read.csv(file = Ed_path, header = TRUE, skip = 9)
  Ed_sim_interp = Hmisc::approxExtrap(x = Ed_sim$Wavelength, y = Ed_sim$Ed_total.W.m.2.nm.,
                         xout = wavelength, method = "linear")$y
  if (dg_comsposite == TRUE) {
    
    abs_cdom_nap_443 = dg_443
    print("Composite absorption for a_dg is provided from user subroutine, i.e. QAA")
    
  } else {
    abs_cdom_nap_443 = abs_cdom_443 + abs_nap_443
    print("Seperate absorption for a_dg is provided from user inputs")
  }
  
  
  Lf_685 <- phi_f * coeff_x[1] * c_chl / (1 + coeff_x[2] * abs_cdom_nap_443 + coeff_x[3] * c_chl)
  
  #Lf_685= 0.001 * Lf_685  #convert the unit 1/??m to 1/nm; Ed: (W/m^-2 nm-1)  v.s. 
                          #Lf:  (W/m^-2 ??m-1) * (mg/m3)-1
 
  # The spectral shape of fluorescence, Lf, was modeled as a Gaussian spectral profile 
  #centered at 685 nm, having a full width at half maximum (FWHM) of 25 nm
  
  Lf_distribute <- function(wv) {
    exp(-4 * log(2) * ((wv - 685)/25)^2) + 0.3 * exp(-4 * log(2) * ((wv - 730)/50)^2)
  }
  
  # Calculate Lf based on the spectral profile and the wavelength
  Lf <- Lf_685 * Lf_distribute(wavelength)
  Rrs_sicf = Lf/(Ed_sim_interp*100)
  return(Rrs_sicf)

}

qaa_output = QAA.v5(waves = wavelength, Rrs = insitu.data)

#Test the SICF function with user supplied values for a_dg
rrs_sicf_dg_input = Rrs_Fluorescence(dg_comsposite=FALSE, dg_443 = qaa_output$a_dg_443,
                            #c_chl = 10,
                            c_chl = Fit.input$chl,
                            abs_cdom_443 = Fit.input$acdom.440,
                            abs_nap_443 = Fit.input$anap.440,phi_f = 0.02)

#Test the SICF function with for a_dg values ibtained from QAAv5
rrs_sicf_dg_qaa = Rrs_Fluorescence(dg_comsposite=TRUE, dg_443 = qaa_output$a_dg_443,
                            #c_chl = 10,  
                            use_analytic_Ed = T,
                            c_chl = Fit.input$chl,
                            abs_cdom_443 = Fit.input$acdom.440,
                            abs_nap_443 = Fit.input$anap.440,phi_f = 0.02)

rrs_forward = Saber_forward(chl = 10, acdom440 = Fit.input$acdom.440, 
                            anap440 =Fit.input$anap.440 , bbp.550 = Fit.input$bbp.550, 
                            verbose = T,z = zB, realdata.exist = FALSE,
                            rb.fraction = fA.set)

rrs.el <- rrs_forward[[1]]$Rrs #Extract AM03 modeled Rrs for elastic component

plot(wavelength, rrs.el, pch=19, col="black")
lines(wavelength, rrs.el+rrs_sicf, col="green")
#lines(wavelength, rrs.el+rrs_sicf, col="orange") #Test effect of quantum yield


#parametric formula to retrieve spectral slope of CDOM + NAP
s_cdm = 0.015 + (0.002/(0.6 + (obsdata[which.min(abs(wavelength - 443))]/obsdata[which.min(abs(wavelength - 555))])))

#parametric formula to retrieve spectral slope of bbp
eta_bbp = 2*(1-(1.2*exp(-0.9 * (obsdata[which.min(abs(wavelength - 443))]/obsdata[which.min(abs(wavelength - 555))]))))


#==================================================================================
##=---------------full model : Light & Water models to be implemented--------------
#==================================================================================

# phi_f <- 0.01  #photons emitted x [photons absorbed]^-1   e.g., 1%.  usually 0.2-1.5%

# Cf <- 1  #the fraction of this emission within 1 nm at 685 nm to the total emission (1/Cf )

# Qa <- 0.5  #fluorescence emitted within the cell with partialy fraction leaving the cell

# a_phy <- list()  #the absorption coefficient of phytoplankton at ??

# #the flux absorbed by phytoplankton
# E0_ <- 1  #(??mol . photon. m-2. s-1. nm-1) is the scalar irradiance 
#(i.e., the excitation irradiance) just below the surface

# K_lambda <- 0.2  #the attenuation of scalar irradiance

# K_Lu <- 0.5 #the attenuation of upwelling irradiance, KLu(685), accounts 
#for the attenuation of upwelling fluorescence radiance

# Lf_685 <- 0.54 * 1/(4*pi) * phi_f/Cf * Qa *  ...
##=---------------full model : Eq. 7.5  , to be implemented------------------------

##test using Eq. 7.23
#Lf_685 <- 0.15 * v_chl.reshape([-1,1])  / (1 +  0.2 * v_chl.reshape([-1,1]) )

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
wrfCDOM <- function(wave1, wave2) {
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
    
    print(paste0("The current ex wavelength = ", wave1, " is not affected by the fDOM wavelength redistribution function"))
    A0 = 0.0
    wrfcdom = 0.0
    
  }

  return(wrfcdom)
}
#----------------------------------------------------------------------------------------
#Discretization to perform the integration for the wavelength redistribution functions
wave.redist.inelastic <- function(wavelength = seq(400,800,10),
                                  include.chl = T,
                                  include.CDOM = T, quant_y_phi = 0.02) {
  
  deltaw <- 1
  waveb <- wavelength
  Nwave <- length(waveb) - 1
  
  if (include.chl == TRUE) {
    
    iChlfl <- 1
  } else {
    
    iChlfl <- 0
  }
  
  if (include.CDOM == TRUE) {
    
    iCDOMfl <- 1
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

#-------------------------------------------------------------------------------------------------------
#Test the wavelength redistribution function with discretization
a = wave.redist.inelastic(quant_y_phi = 0.05) #call this variable later to do the final integral

plot(wavelength[-1], a[["sicf"]][1,], pch=19, type="l", lwd=2.5, col="green",
     ylab = expression(paste("Fluoroscence Efficieny (1/nm)")), xlab = "Wavelength" )
lines(wavelength[-1], a[["fDOM"]][1,], col="goldenrod2", lwd=2.5) #Show 

legend(x = "topleft",          # Position
       legend = c("SICF","fDOM"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = rev(c("goldenrod2", "green")),           # Line colors
       lwd = 2)

#-------------------------------------------------------------------------------------------------------
#Implement the final source scattering function for inelastic scattering 

######Scalar irradiance calculation#####

#---------------------------------------------------------------------------
#Get Ed0 from HL simulated Ed data
#---------------------------------------------------------------------------
Ed_path = "./input-spectra/Ed_HL.csv"
Ed_sim = read.csv(file = Ed_path, header = TRUE, skip = 9)


Ed_sim_interp = Hmisc::approxExtrap(x = Ed_sim$Wavelength, y = log(Ed_sim$Ed_total.W.m.2.nm.),
                                    xout = wavelength, method = "linear")$y
Ed_sim_interp = exp(Ed_sim_interp)*100 #extrapolate and change units

sun_view = snell_law(sun = sun, view = view) #calculate sun zenith

E0_sim_interp = Ed_sim_interp/(2*cos(sun_view$sun_w)) #Convert to scalar irradiance

E0_sim_interp = 0.96*E0_sim_interp #Convert to underwater E0 (Check with Simon)

#---------------------------------------------------------------------------
#Get the Ed from COP-S package calculated value
#---------------------------------------------------------------------------
load("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/COPS_Kildir/BIN/WISE_CAST_001_190818_205712_URC.csv.RData")

Kd =cops$K.EdZ.surf; #Kd for Ed at surface

Edz_0m_loess = cops[["EdZ.0m"]] #Ed_0- 
Edz_0m_linear = cops$EdZ.0m.linear #(same with above, the extrapolation is different)

Ed0_0p = cops$Ed0.0p #Ed_0+

ed_waves = cops$EdZ.waves #Ed waves for COP-S

#Do the interpolation to simulation wavelengths
kd_interp = exp(approx(x = ed_waves, y= log(Kd), xout = wavelength, method = "linear")$y)

Edz_0m_loess_interp = exp(approx(x = ed_waves, y= log(Edz_0m_loess), 
                                 xout = wavelength, method = "linear")$y)

Edz_0m_linear_interp = exp(approx(x = ed_waves, y= log(Edz_0m_linear), 
                                  xout = wavelength, method = "linear")$y)

Ed0_0p_interp = exp(approx(x = ed_waves, y= log(Ed0_0p), 
                           xout = wavelength, method = "linear")$y)

#---------------------------------------------------------------------------
#Calculate E0 from absorption and Kd from in situ data following Morel 1991
#---------------------------------------------------------------------------
iop_data = read.surface.IOPs.wise(station.args = "OUT-F18") #load station specific IOP and AOP
a_tw = data.frame(test_IOP_func$a_data) #extract non-water absorption
a_tw$at_w[a_tw$at_w <= 0] <- 1e-10

a_tw_interp = exp(approx(x=as.numeric(as.character(a_tw$wave)), 
                         y = log(as.numeric(a_tw$at_w)), xout = wavelength, 
                         method = "linear")$y) #interpolate the non-water absorption

a_tw_interp[is.na(a_tw_interp)] = 1e-10 #replace -ve values with very small ones

E0_morel <- (kd_interp* Ed_interp)/(a_tw_interp)

# E0_morel #Morel given E0
# E0_interp #COp-S E0 (not SURE)
# Ed0_m_interp #COp-S E0 (not SURE) 

#---------------------------------------------------------------------------
#Calculate the Ed using Gergg & Carder model
#---------------------------------------------------------------------------
GreggCarder.data()

require(lubridate)
#x = as.Date('2019-08-18')
x = as.Date(cops$date.mean)
jday_no = yday(x)

test_Ed = GreggCarder.f(the = cops$sunzen, 
                        lam.sel = seq(400,800,10),
              jday = jday_no, rlon = cops$longitude, rlat = cops$latitude, debug = T)

Ed0.0p = test_Ed$Ed
Ed0_dir.0p = test_Ed$Edir
Ed0_dif.0p = test_Ed$Edif

plot(wavelength, Ed0_0p_interp, pch=19, ylim=c(0,150))
lines(wavelength, Ed0.0p*100, col="black", lwd=2)


rhoF <- GreggCarder.sfcrfl(rad = 180.0/pi, theta=cops$sunzen, ws=5)
sun_view_op = snell_law(view = 0 , sun = cops$sunzen)

fEd_dir = 0.8
#=================================================================================
rh = 60

# hour_angle = 15 * (local_time - 12)
# 
# 
# alt = asin(sin(lat) * sin(zen) + cos(lat) * cos(zen) * cos(hour_angle))
# 
# 
# fEd_dir <- function(V, RH, theta, alpha) {
#   lnV <- log(V)
#   lntheta <- log(theta)
#   lnalpha <- log(alpha)
#   return(0.356 - 0.128*lnV + 0.044*RH - 0.245*lntheta - 0.777*lnalpha)
# }
# 
# # Example usage
# fEd_dir(V = 15, RH = rh, theta = 42, alpha = 0.2)


#fEd_dir <- 1 - 1.04 * exp(-1.45/cos(sun_view_op$sun_w)) * (1 - exp(-0.63/cos(sun_view_op$sun_w))) * (1 - exp(-0.41*RH/100))

# Print result
print(fEd_dir)


#================================================================================


# 
# Ed0.0m = (Ed0.0p * fEd_dir       * (1 - rhoF$rod)) + # direct
#          (Ed0.0p * (1 - fEd_dir) * (1 - rhoF$ros))   # diffuse

Ed0.0m = (Ed0_dir.0p *(1 - rhoF$rod)) + # direct
        (Ed0_dif.0p * (1 - rhoF$ros))   # diffuse

lines(wavelength, Edz_0m_linear_interp, col="green", lwd=2)
lines(wavelength, Ed0.0m*100, col="green", lwd=2, lty="dashed")

#Convert it to scalar irradiance
sun_view_om = snell_law(view = 0 , sun = sunzen_below(cops$sunzen))
Ed0m_0 = (Ed0.0m*100) / cos(sun_view_om$sun_w)

lines(wavelength, Ed0m_0, col="goldenrod2", lwd=2)

#######Calculate the CDOM absorption which is the scattering function for fDOM#####
abs_CDM <- rep(0,length(wavelength))
abs_CDM_440 = Fit.input$acdom.440
S_CDM = 0.014

for (i in 1:length(wavelength)){
  
  abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(wavelength[i] - 440))
  
}


#######Do the integration with respect to wavelength of the final source function of fDOM#####
# the equation to be integrated is:a_CDOM*f_F(ex, em)
waveb <- seq(400,800,10)
Nwave <- length(waveb) - 1

WRF_CDOM = a$fDOM

sum0 = matrix(0, ncol = Nwave, nrow = Nwave)

for (jwave in 2:Nwave) {
  
  for (iwave in 1:(jwave-1)) {
    
    #Eozi <- sum(E0_sim_interp[1:iwave]) #Get the scalar irradiance at iwave
    #Eozi <- E0_sim_interp[iwave]
    
    #Eozi <- E0_morel[iwave]
    
    #Eozi <- E0_interp[iwave]
    
    Eozi = sum(Ed0m_0[1:iwave])
    
    
    absCDOM <- sum(abs_CDM[1:iwave]) #Get the CDOM absorption at iwave
    #absCDOM <- abs_CDM[iwave]
    
    sum0[iwave,jwave] <- sum0[iwave,jwave] + (Eozi * absCDOM * WRF_CDOM[iwave, jwave])
    
  }
}

fDOM_rad = sum0[1,]/(4*pi) #add the isotropic distributed phase function of fDOM (*(1/4*pi))

plot(wavelength[-length(wavelength)], fDOM_rad, pch=19, col="goldenrod2")

# fDOM_rrs_ed_sim = fDOM_rad / Ed_sim_interp[-length(Ed_sim_interp)]
# 
# fDOM_rrs_ed_cops = fDOM_rad / Ed_interp[-length(Ed_interp)]
# 
# fDOM_rrs_ed_morel = fDOM_rad / Ed0_m_interp[-length(Ed0_m_interp)]

fDOM_rrs_ed_analytical = fDOM_rad / (Ed0.0m[-length(Ed0.0m)]*100)

####Plot the fDOM equivalent rrs#################
plot(wavelength, surface_rrs_translate(Rrs = insitu.data), ylim=c(0, 0.005),
     pch=19, col="black", ylab = expression(paste(R["rs"], "(",0^"-",")")))

rrs_final = rrs.forward.am.param.conc.true_iop_sicf[-length(rrs.forward.am.param.conc.true_iop_sicf)] + fDOM_rrs_ed_analytical

output_rrs = rrs.forward.am.param.conc.true_iop_sicf

output_rrs[-length(output_rrs)] = rrs_final

lines(wavelength, output_rrs,
      col="green", lwd=2)
lines(wavelength, rrs.forward.am.param.conc.true_iop , col="navyblue", lwd=2)

legend(x = "topright",          # Position
       legend = c("insitu","elastic+ inelastic","elastic"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("black","green", "blue"),           # Line colors
       lwd = 2, bty = "n")

#-------------------------------------------------------------------------------------------
# for (jwave in 2:Nwave) {
#   
#   for (iwave in 1:(jwave-1)) {
#     
#     #sum = sum + (WRF_CDOM[iwave,jwave] * abs_CDOM * ((kd_interp * Ed_interp)/abs_CDOM))
#     sum = sum + (WRF_CDOM[iwave,jwave] * abs_CDOM * (E0_interp))
#     
#   }
# }

#-------------------------------------------------------------------------------------------
# library(pracma)
# 
# deltaw <- 1
# waveb <- seq(400,800,10)
# Nwave <- length(waveb) - 1
# 
# WRF_Chl <- matrix(0, nrow = Nwave, ncol = Nwave)
# 
# for (j in 2:Nwave) {
#   wj <- waveb[j] + 0.5 * deltaw
#   nj <- round((waveb[j + 1] - waveb[j]) / deltaw)
#   factor <- deltaw * deltaw / (waveb[j + 1] - waveb[j])
#   
#   for (i in 1:(j-1)) {
#     wi <- waveb[i] + 0.5 * deltaw
#     ni <- round((waveb[i + 1] - waveb[i]) / deltaw)
#     
#     sumchl <- 0
#     for (jj in 1:nj) {
#       wavej <- wj + (jj - 1) * deltaw
#       for (ii in 1:ni) {
#         wavei <- wi + (ii - 1) * deltaw
#         sumchl <- sumchl + wrfChl(wavei, wavej)
#       }
#     }
#     
#     WRF_Chl[i, j] <- factor * sumchl
#   }
# }

# WRF_Chl <-matrix(0, nrow = Nwave, ncol = Nwave)
# WRF_CDOM <- matrix(0, nrow = Nwave, ncol = Nwave)
# WRF_Raman <- matrix(0, nrow = Nwave, ncol = Nwave)
# 
# for (j in 2:Nwave) {
#   wj <- waveb[j] + 0.5 * deltaw
#   nj <- as.integer((waveb[j + 1] - waveb[j])/deltaw)
#   factor <- deltaw * deltaw / (waveb[j + 1] - waveb[j])
#   
#   for (i in 1:(j - 1)) {
#     wi <- waveb[i] + 0.5 * deltaw
#     ni <- as.integer((waveb[i + 1] - waveb[i])/deltaw)
#     
#     if (iChlfl != 0) {
#       sumchl <- 0
#       for (jj in 1:nj) {
#         wavej <- wj + as.numeric(jj - 1) * deltaw
#         for (ii in 1:ni) {
#           wavei <- wi + as.numeric(ii - 1) * deltaw
#           sumchl <- sumchl + wrfChl(wavei, wavej)
#         }
#       }
#       
#       WRF_Chl[i, j] <- factor * sumchl
#     }
#     
#     if (iCDOMfl != 0) {
#       sumcdom <- 0
#       for (jj in 1:nj) {
#         wavej <- wj + as.numeric(jj - 1) * deltaw
#         for (ii in 1:ni) {
#           wavei <- wi + as.numeric(ii - 1) * deltaw
#           sumcdom <- sumcdom + wrfCDOM(wavei, wavej)
#         }
#       }
#       
#       WRF_CDOM[i, j] <- factor * sumcdom
#     }
#   }
# }
#-----------------------------------------------------------------------------------------------
