#==============================================================================================
# main.R performs a series of operations on Radiative transfer in optically complex case-II
# waters. Some of the major functionalities are listed below:

# 1. Model the Remote-Sensing Reflectance in optically complex waters from user-given [chl],
# CDOM and NAP absorption at 440nm wavelength and in case of shallow water, depth and 
# aerial fraction of bottom types using SA models. Two types of SA models; I.Albert & Mobley 03,
# II. Lee et al. 98. The AM03 model has two modes for the calculation of slopes in BGC models,
# i.e.constant and dynamically obtained from rrs. The AM03 model is also able to model SICF if
# quantum yield is provided.

# 2. Retrieve [chl], a_cdom+nap(440) and z with Rb (shallow) from user given rrs by inverting
# SA model. User can opt between the AM03 and Lee98 SA models to invert. The inversion can be
# performed by finite difference based optimization, constrained gradient based optimization,
# Stochastic optimization and Markov Chain Monte Carlo optimization. For gradient based method,
# an analytical jacobian can be supplied to drastically improve the optimization speed.

# Note, all methods except the Stochastic Optimization and MCMC are sensitive to initial values.
# A prefit of parameters is available following the CRISTAL LUT approach of Mobley 200~. For
# gradient based methods, it is advised to perform the prefit prior to optimization.

#3. The AM03 model for shallow waters has two mode of inversion, I. unconstrained: all parameters
# are unknown, II. constrained : only depth and bottom albedo are unknown. for the bottom
# albedo, there are two mode, I. retrieve as a weighted sum of fractions of known type of n classes,
# II. Retrieve the entire spectral Rb (only possible if all other parameters are known).
 

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#==============================================================================================
rm(list=ls(all=TRUE))
setwd("/home/musk0001/R_inverse_wasi")

source("./saber_forward.R")
source("./saber_forward_grad.R")
source("./saber_forward_paramteric.R")
source("./saber_forward_parametric_conc.R")
source("./saber_forward_parametric_conc_wise.R")
source("./SABER.sicf.R")
source("./QAA_v5.R")
source("./saber_inverse.R")
source("./lee_forward.R")

source("./snell_law.R")
#source("C:/R/Cops/R/GreggCarder.R")

source("./retrive.rb.spectral.wise.R")
source("./generate.rb.spectral.R")
source("./read.surface.IOPs.wise.R")
source("./inelastic.wavelength.redistribution.R")

source("./solve.objective.inverse.R")
source("./solve.objective.inverse.shallow.R")
source("./solve.objective.inverse.shallowv2.R")
source("./solve.objective.inverse.shallow.constrained.batch.R")
source("./mcmc.functions.R")
#-------------------------------------------------------------------------------------------
library(dplyr)
library(readxl)
library(stats4)
library(MASS)
library(dglm)
library(fitdistrplus)
library(Riops)
library(Cops)
library(ggplot2)
library(rho)
library(marqLevAlg)
library(BayesianTools)
library(coda)
library(stats)
library(Deriv)
library(mcmc)
library(plot3D)
library(limSolve)
library(nloptr)
library(minpack.lm)
library(pracma)
#-------------------------------------------------------------------------------------------
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

sunzen_below = function(sun_zen_aove=45){
  sun_zen_below_rad = asin(sin(sun_zen_aove*(pi/180))/1.333)
  sun_zen_below_deg = sun_zen_below_rad*(180/pi)
  return(sun_zen_below_deg)
}
#--------------------------------------------------------------------------
## 1. Inputs
#--------------------------------------------------------------------------

#1.1 Water type specifications
type_case_water = 2
type_Rrs_below = "deep"
type_Rrs_below_rb = "shallow"
type_Rrs_water = "below_surface"

#1.2 Data and station specifications
batch=FALSE #Set TRUE for IOCCG dataset duplication; Set FALSE for user wanted inputs

insitu.present=TRUE #Set TRUE if actual in situ or simulated observations exist; 
                     #else set FALSE
statname = "MAN-R01" #if insitu.present = TRUE, set the station-name 

insitu_type = c("HL", "COPS") 
insitu.type = insitu_type[2] #<<USER INPUT >> selection for type of in situ data

HLpath_el = "/home/musk0001/L2/20190818_StationOUT-F18/IOP/hl_simulation/HL_el_table_T-F18.csv" #<<USER PATH>>

HLpath_inel = "/home/musk0001/L2/20190818_StationOUT-F18/IOP/hl_simulation/HL_inel_table_T-F18.csv" #<<USER PATH>>

use_bb_nup = TRUE #Set TRUE If bbp555 should be calculated using in situ spectral slope 

if (batch == TRUE) {
  j=376 #  <<USER INPUT >> (No. of IOCCG data point among 500) (for single run mode)
        #For batch inversion of all data in the IOCCG dataset refer to the R code 
        #"Saber.batch.deep.R" and "Saber.batch.deep.userdefined.R" (user defined randomly 
        #splits the data into user defined partitions and inverts them respectively.)
}


plot=FALSE #Set TRUE when the ggplot2 outputs needed to be saved onto disk


#1.3 Inversion specifications
preFit = FALSE

use.lklhood= TRUE
use.wise.prior = FALSE 
use.nomad.prior = TRUE
use.ioccg.prior = FALSE

#1.4  Desired Wavelength for the simulation
wavelength <- seq(400,800,10)

#=======================================================================================================
#--------------------------DATA HANDLING MANUAL: START----------------------------------
# Observed/in situ Rrs for following parameters to be input in RT
# Rrs.cops <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Cops_table_T-F18.csv", 
#                     header = T)
# Rrs.albert <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Albert_table_T-F18.csv",
#                        header = T)
# Rrs.HL <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/HL_el_table_T-F18.csv",
#                        header = T)
# 
# Rrs_obs <- as.numeric(Rrs.albert[1,-1])
# Rrs_obs_wl <- c(320, 330, 340, 380, 412, 443, 465, 490, 510, 
#                 532, 555, 589, 625, 665, 683, 694, 710, 780, 875)
# 
# Rrs_obs.interp = Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs, xout = wavelength,
#                                      method = "linear")$y
# 
# Rrs_obs.cops <- as.numeric(Rrs.cops[1,-1])
# Rrs_obs.cops.interp = Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs.cops, xout = wavelength,
#                                      method = "linear")$y
#--------------------------DATA HANDLING MANUAL: END----------------------------------
#=======================================================================================================

#1.5 Viewing geometry in degrees
view = 0
sun_above  = 60
sun = sunzen_below(sun_zen_aove = sun_above) #Make sure you pass "sun" for subsurface calculation

#1.6 bottom depth
zB=2

#1.7 Areal fraction of bottom reflectance
fA0=0; # constant (Excluded from Model) 
fA1=0.2; # sand
fA2=0.2; # sediment
fA3=0; # Chara contraria
fA4=0.6; # Potamogeton perfoliatus
fA5=0; # Potamogeton pectinatus

#fA.set= c(fA0,fA1,fA2,fA3,fA4,fA5)
fA.set= c(fA1,fA2,fA3,fA4,fA5) #exclude the constant fA0

source("./retrive.rb.spectral.R")
source("./generate.rb.spectral.R")

#1.7 Atmospheric conditions

# Irradiance intensities [1/sr]
g_dd=0.05; g_dsr=0; g_dsa=0;

# Intensities of light sources 
f_dd= 1; f_ds= 1;

# Angstrom exponent
alpha = 1.317;

# Atmospheric pressure 
P = 1013.25; # [mbar]

# Relative Humidity
RH = 0.60;

# Scale height for ozone
Hoz=0.300; # [cm]

# Scale height of the precipitate water in the atmosphere
WV= 2.500; # [cm]

#-----------------------------------------------------------------------------------------------------
#2. Instantiate IOPs and AOPs into working environment 
#-----------------------------------------------------------------------------------------------------

##2.1  Read the IOCCG database 

#(Note, this dataset is totally driven by [chl]; a_phi, a_d, a_g and bbp, all 
#are derived using empirical formulations with [chl], hence, individual a and bbp 
#values may yield different rrs from the dataset)

#To see the difference between [chl] driven parameter based simulation and free IOP parameter
#based simulation refer to the following code till prior to Section 2.2. 

if (batch == TRUE) {
  insitu.data.HL <- read_excel_allsheets(paste0(getwd(),"/IOP_AOP_Sun60.xls"))
  #water IOP
  water.data <- insitu.data.HL$Basics[6:46, 1:3]                           
  names(water.data) <- c("wavelength", "a_w", "bb_w") 
  
  a.ph = insitu.data.HL$a_ph; a.g = insitu.data.HL$a_g; a.d = insitu.data.HL$a_dm
  bb = insitu.data.HL$bb_ch + insitu.data.HL$bb_dm
  
  #chl, acdom440, anap440 & bbp550 from IOCCG dataset 
  chldata <- insitu.data.HL$a_ph[-1,1]
  
  acdom440data <- insitu.data.HL$a_g[which(names(insitu.data.HL$a_g[1,]) == "440")]
  
  anap440data <- insitu.data.HL$a_dm[which(names(insitu.data.HL$a_dm[1,]) == "440")]
  
  bbp550data <- insitu.data.HL$bb[which(names(insitu.data.HL$bb[1,]) == "550")] -
    as.numeric(water.data$bb_w[water.data$wavelength == 550])
  
  HL.deep.iop <- data.frame("chl"=chldata, "acdom440"= acdom440data$`440`,
                            "anap440"= anap440data$`440`, "bbp550"= bbp550data$`550`)
  #sub-surface (0^-) Rrs
  rrs.HL <- as.matrix(insitu.data.HL$r_rs); rrs.HL.wl <- as.numeric(names(insitu.data.HL$r_rs))
  
  #Generate the AM03 simulated Rrs(0^-) using loaded IOPs
  rrs.forward.SABER <- matrix(nrow = length(HL.deep.iop$chl), ncol = length(wavelength),0)
  
  #Create the Progress Bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(HL.deep.iop$chl), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")
  
  
  for (i in 1:length(HL.deep.iop$chl)) { #Create forward SABER LUT
    temp1 <- as.numeric(HL.deep.iop[i,])
    temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                           anap440 =temp1[3], bbp.550 = temp1[4], verbose = F, 
                           realdata = rrs.HL[i,])
    rrs.forward.SABER[i,] <- temp2[[1]]$Rrs
    setTxtProgressBar(pb, i)
    if (i == length(HL.deep.iop$chl)) {
      cat(paste0("\033[0;32m","#LUT generation for AM03 FINISHED#","\033[0m","\n"))
    }
    #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.SABER) - i), " remaining","\033[0m","\n"))
    
  }
  
  #Generate the LEE99 simulated Rrs(0^-) using loaded IOPs
  rrs.forward.lee <- matrix(nrow = length(HL.deep.iop$chl), ncol = length(wavelength),0)
  
  for (i in 1:length(HL.deep.iop$chl)) { #Create forward SABER LUT
    temp1 <- as.numeric(HL.deep.iop[i,])
    temp2 <- Lee_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = temp1[4], verbose = F, realdata = rrs.HL[i,] )
    rrs.forward.lee[i,] <- temp2[[1]]$Rrs
    #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.lee) - i), " remaining","\033[0m","\n"))
    setTxtProgressBar(pb, i)
    if (i == length(HL.deep.iop$chl)) {
      cat(paste0("\033[0;32m","#LUT generation for LEE99 FINISHED#","\033[0m","\n"))
    }
  }
  
  #Test the 1:1 scatter b/w SABER forward and HL
  plot(as.matrix(rrs.HL), as.numeric(rrs.forward.lee), xlim=c(0,0.05), ylim=c(0,0.05),
       xlab="Hydrolight", ylab="Modelled")
  abline(0,1, col="green", lwd=3)
}

#------------------------------------------------------------
#2.2 Obtain in vivo values for IOP and BGC from sampling data
#------------------------------------------------------------
if (insitu.present == TRUE & !is.null(insitu.type)) {
  invivo_data = get_in_situ_params(station_name = statname, use_bb_nup = TRUE)
}

#---------------------------------------------------------
#3. prepare the relevant data structure
#---------------------------------------------------------
if (insitu.present == TRUE & batch == TRUE) { #Set the observed Rrs manually from IOCCG
  
  insitu.data <-rrs.HL[j,]
} 

if (insitu.present == TRUE & batch == FALSE & insitu.type == "COPS"){#Set the observed Rrs
                                                                  #manually from field data
  
  IOP_AOP_surf = suppressWarnings(read.surface.IOPs.wise(station.args = statname))
  insitu.data = IOP_AOP_surf$Rrs_0p
} 

if (insitu.present == TRUE & batch == FALSE & insitu.type == "HL") {#Set the observed Rrs 
                                                #manually from HL simulation of field data
  
  Rrs.hl <- read.csv(HLpath_inel, header = T)
  
  Rrs_obs <- as.numeric(Rrs.hl[1,-1])
  Rrs_obs_wl <- c(320, 330, 340, 380, 412, 443, 465, 490, 510,
                  532, 555, 589, 625, 665, 683, 694, 710, 780, 875)
  
  Rrs_obs.interp_HL = approx(Rrs_obs_wl, Rrs_obs, xout = wavelength,
                          method = "linear")$y
  insitu.data = Rrs_obs.interp_HL
} 

# plot(wavelength, surface_rrs_translate(insitu.data))
# lines(wavelength, Rrs_obs.interp_HL, lwd=2, col="red")
# lines(wavelength, rrs.forward.am.param.conc.true_iop, lwd=2, col="navyblue")
# 
# lines(wavelength, Rrs_obs.interp_HL, lwd=2, col="red", lty="dashed")
# lines(wavelength, rrs.forward.am.param.conc.true_iop_sicf_fdom, lwd=2, col="navyblue", lty="dashed")
# 
# plot(wavelength, surface_rrs_translate(insitu.data))
# nc <- 2   
# Vars <- c("HL", "SA model")
# nr <- length(Vars)
# 
# legend("bottom", rep(Vars, nc), col = c("red", "navyblue"), 
#        lty = rep(1:nc, each = nr), ncol = nc, cex = 0.8, 
#        title = "Elastic              Elastic+inelastic", bty = "n")


if (insitu.present == FALSE & batch == FALSE) {#Set the observed Rrs as a QSAA model retrieved Rrs
  
  rrs.demo <- read.csv("./input-spectra/demo_rrs.csv",header = T)
  
  
  rrs.demo.Om <- as.numeric(rrs.demo[1,-1])
  rrs.demo.wl <- c(320, 330, 340, 380, 412, 443, 465, 490, 510,
                  532, 555, 589, 625, 665, 683, 694, 710, 780, 875)
  
  rrs.demo.Om.demo = approx(rrs.demo.wl, rrs.demo.Om, xout = wavelength,
                          method = "linear")$y
  insitu.data = rrs.demo.Om.demo #Note, this is actually a trick to let the forward model run
                                #however, this is not the right insitu.data 
} 

#-----------------------------------------------------------
##4. Single FORWARD RUN
#-----------------------------------------------------------

#4.1 Create input vector
if (batch == TRUE) { 
  #For IOCCG
  Fit.input <- data.frame("chl"=  HL.deep.iop$chl[j],
                          "acdom.440"=  HL.deep.iop$acdom440[j],
                          "anap.440"= HL.deep.iop$anap440[j],
                          "bbp.550" = HL.deep.iop$bbp550[j])
} 

if (insitu.present == TRUE) {
  #Data from in-vivo vals
  Fit.input <- data.frame("chl"=  invivo_data$chl_invivo, 
                          "acdom.440"=  invivo_data$abs_invivo$a_cdom ,
                          "anap.440"= invivo_data$abs_invivo$a_nap ,
                          "bbp.550"=invivo_data$bbp555)
  
} 

if (insitu.present == FALSE) {
  #Manual entries
  Fit.input <- data.frame("chl"=  4.96, #Data for OUT-F18 from in vivo vals
                          "acdom.440"=  0.9322 ,
                          "anap.440"= 0.07 ,
                          "bbp.550"=0.0072)
  
}

#4.2 Auxiliary declarations
base.bbp = Fit.input$bbp.550
type_Rrs_below_jacobian = "deep"
Cops::GreggCarder.data()

IOP_files = list.files("./Rb_spectral/surface_iops/", full.names = T)

idx_a = grep(IOP_files, pattern = paste0("abs_surf_",statname))
idx_bb = grep(IOP_files, pattern = paste0("bb_surf_",statname))


#-----------------------------------------------------------------------------
#4.3 Test different bio-optical parametrization effect on rrs in forward SABER
#-----------------------------------------------------------------------------

#4.3.1. Constant values for s_g and s_nap (a_d,g spectral slopes) and \eta (bbp spectral slope)
#with no spectral shape normalization and no SICF
forward.op.am <- Saber_forward(chl = Fit.input$chl, acdom440 = Fit.input$acdom.440, 
                            anap440 =Fit.input$anap.440 , bbp.550 = Fit.input$bbp.550, 
                            #realdata = insitu.data, 
                            z = zB, rb.fraction = fA.set,
                            verbose = F, realdata.exist = TRUE,
                            realdata = surface_rrs_translate(Rrs = insitu.data),
                            plot = F)

rrs.forward.am <- forward.op.am[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.2. parametric values for s_dg (a_dg slope) and \eta (bbp slope) following QAAv5
#with no spectral shape normalization and no SICF
forward.op.am.param <- Saber_forward_paramteric(chl = Fit.input$chl, 
                           acdom440 = Fit.input$acdom.440, 
                           anap440 =Fit.input$anap.440 , 
                           bbp.550 = Fit.input$bbp.550,
                           #realdata = rrs.forward.am,
                           realdata.exist = T,
                           realdata = surface_rrs_translate(Rrs = insitu.data),
                           dg_composite = T,
                           slope.parametric = T,
                           verbose = F,
                           sicf = F, plot = F)


rrs.forward.am.param <- forward.op.am.param[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.3.1. parameteric values for s_dg (a_dg slope) and \eta (bbp slope) following QAAv5
#with spectral shape normalization and no SICF
forward.op.am.param.conc.dg_sep <- Saber_forward_paramteric_conc(chl = Fit.input$chl, 
                                   acdom440 = Fit.input$acdom.440, 
                                   anap440 =Fit.input$anap.440 ,
                                   bbp.550 = Fit.input$bbp.550, 
                                   a_dg = NULL,
                                   #realdata = rrs.forward.am,
                                   realdata = surface_rrs_translate(Rrs = insitu.data),
                                   slope.parametric = F,
                                   dg_composite = F,
                                   use_spectral_shape_chl = F,
                                   use_spectral_shape_dg = T,
                                   sicf = F, q_phi = 0.02, 
                                   verbose = F, plot = F)

rrs.forward.am.param.conc.dg_sep <- forward.op.am.param.conc.dg_sep[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.3.2. parametric values for s_dg (a_dg slope) and \eta (bbp slope) following QAAv5
#with spectral shape normalization and no SICF
forward.op.am.param.conc.dg_comp <- Saber_forward_paramteric_conc(chl = Fit.input$chl, 
                                    acdom440 =NULL, 
                                    anap440 =NULL , 
                                    a_dg = Fit.input$acdom.440 + Fit.input$anap.440,
                                    bbp.550 = Fit.input$bbp.550,
                                    #realdata = rrs.forward.am,
                                    realdata = surface_rrs_translate(Rrs = insitu.data),
                                    slope.parametric = T,
                                    dg_composite = T,
                                    use_spectral_shape_chl = F,
                                    use_spectral_shape_dg = T,
                                    sicf = F, q_phi = 0.02, 
                                    verbose = F, plot = F)

rrs.forward.am.param.conc.dg_comp <- forward.op.am.param.conc.dg_comp[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.3.3. parametric values for s_dg (a_dg slope) and \eta (bbp slope) following QAAv5
#with spectral shape normalization and SICF+fDOM
forward.op.am.param.conc.dg_comp_sicf_fdom <- Saber_forward_paramteric_conc_wise(
                                              use_true_IOPs = F,
                                              
                                              chl = Fit.input$chl, 
                                              acdom440 = NULL, 
                                              anap440 = NULL, 
                                              a_dg =  Fit.input$acdom.440 + Fit.input$anap.440 ,
                                              bbp.550 = Fit.input$bbp.550,
                                              
                                              #realdata = rrs.forward.am,
                                              realdata = surface_rrs_translate(Rrs = insitu.data),
                                              
                                              slope.parametric = T,
                                              dg_composite = T,
                                              use_spectral_shape_chl = F,
                                              use_spectral_shape_dg = T,
                                              
                                              sicf = T, q_phi = 0.05, 
                                              
                                              fDOM = T,
                                              sunzen_Ed = -99, 
                                              lat_Ed = 49.02487, lon_Ed = -68.37059,
                                              date_time_Ed = "2019-08-18 20:59 GMT", 
                                              Ed_fDOM_path = "./input-spectra/Ed_HL.csv",
                                              use_fDOM_rad = F,
                                              
                                              verbose = F, plot = T)

rrs.forward.am.param.conc.dg_comp_sicf_fdom <- forward.op.am.param.conc.dg_comp_sicf_fdom[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.4.1. full spectral IOPs are provided and no SICF+fDOM
forward.op.am.param.conc.true_iop <- Saber_forward_paramteric_conc_wise(use_true_IOPs = T, 
                                      a_non_water_path = IOP_files[idx_a],
                                      bb_non_water_path = IOP_files[idx_bb],
                                      chl = Fit.input$chl, 
                                      acdom440 =NULL, 
                                      anap440 =NULL , 
                                      a_dg = Fit.input$acdom.440 + Fit.input$anap.440,
                                      bbp.550 = Fit.input$bbp.550,
                                      #realdata = rrs.forward.am,
                                      realdata = surface_rrs_translate(Rrs = insitu.data),
                                      slope.parametric = T,
                                      dg_composite = T,
                                      use_spectral_shape_chl = F,
                                      use_spectral_shape_dg = T,
                                      sicf = F, q_phi = 0.02, 
                                      fDOM = F,
                                      verbose = F, plot = F)

rrs.forward.am.param.conc.true_iop <- forward.op.am.param.conc.true_iop[[1]]$Rrs #Extract AM03 modeled Rrs

#4.3.4.2. full spectral IOPs are provided and included SICF+fDOM
forward.op.am.param.conc.true_iop_sicf_fDOM <- Saber_forward_paramteric_conc_wise(
                                              use_true_IOPs = T, 
                                              a_non_water_path = IOP_files[idx_a],
                                              bb_non_water_path = IOP_files[idx_bb],
                                              
                                              chl = Fit.input$chl, 
                                              acdom440 =NULL, 
                                              anap440 =NULL , 
                                              a_dg = Fit.input$acdom.440 + Fit.input$anap.440,
                                              bbp.550 = Fit.input$bbp.550,
                                              
                                              #realdata = rrs.forward.am,
                                              realdata = surface_rrs_translate(Rrs = insitu.data),
                                            
                                              slope.parametric = T,
                                              dg_composite = T,
                                              use_spectral_shape_chl = F,
                                              use_spectral_shape_dg = T,
                                              
                                              sicf = T, q_phi = 0.05,
                                              
                                              use_analytic_Ed = T,
                                              
                                              fDOM = T, sunzen_Ed = -99,
                                              lat_Ed = 49.02487, lon_Ed = -68.37059,
                                              date_time_Ed = "2019-08-18 20:59 GMT", 
                                              Ed_fDOM_path = "./input-spectra/Ed_HL.csv",
                                              use_fDOM_rad = F,
                                            
                                            verbose = F, plot = T)

rrs.forward.am.param.conc.true_iop_sicf_fdom <- forward.op.am.param.conc.true_iop_sicf_fDOM[[1]]$Rrs #Extract AM03 modeled Rrs


#4.3.5. Plot the Rrs from different Bio-optical parametrizations
forward_rrs_collection = data.frame("wave" = wavelength, "insitu" = surface_rrs_translate(insitu.data),
                                    "saber_param_dg_sep" = rrs.forward.am.param.conc.dg_sep,
                                    "saber_param_dg_comp" = rrs.forward.am.param.conc.dg_comp,
                                    "saber_param_dg_comp_inel" = rrs.forward.am.param.conc.dg_comp_sicf_fdom,
                                    "saber_param_iop" = rrs.forward.am.param.conc.true_iop,
                                    "saber_param_iop_inel" = rrs.forward.am.param.conc.true_iop_sicf_fdom)

xmin = 400; xmax= 800; xstp=100
ymin= 0; ymax=signif(1.25*max(forward_rrs_collection$insitu), digits = 1);ystp= ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(data = forward_rrs_collection)  + 
  
  geom_point(aes(y=insitu,color="xx1", x = wave),size=2,show.legend = F) +
  
  geom_line(aes(y=saber_param_dg_sep, x = wave,color="xx2", linetype="aa1"),
            size=1.3,show.legend = F)+
  geom_line(aes(y=saber_param_dg_comp, x = wave,color="xx3", linetype="aa1"),
            size=1.3,show.legend = F)+
  geom_line(aes(y=saber_param_dg_comp_inel,x = wave,color="xx3", linetype="aa2"), 
            size=1.3,show.legend = F)+
  
  geom_line(aes(y=saber_param_iop, x = wave,color="xx4", linetype="aa1"), 
            size=1.3, show.legend = T)+
  geom_line(aes(y=saber_param_iop_inel, x = wave,color="xx4", linetype="aa2"), 
            size=1.3,show.legend = T)+
  
  scale_colour_manual(labels = c(expression(paste(italic("observed"))),
                                 expression(paste("abs-",italic("a")["dg"],italic(",b")["bp"])),
                                 expression(paste("QAA-",italic("a")["dg"],italic(",b")["bp"])),
                                 expression(paste(italic("field-"),italic("a")["t-w"],italic(",b")["bp"]))),
                      values = c("black","green", "red", "purple")) +
  scale_linetype_manual(labels = c("elastic", "elastic+inelastic"),
                        values = c(1,4))+
  #ggtitle(paste0(stationlist[i])) +
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 15, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.55, 0.98),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.5,0.0,0.0), "cm"),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g  
if (plot == TRUE) {
  ggsave(paste0("./forward.SABER.collection_",statname,".png"), plot = g,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
}


#4.4 Test the Lee SA model
forward.op.lee <- Lee_forward(chl = Fit.input$chl, acdom440 = Fit.input$acdom.440, 
                              anap440 =Fit.input$anap.440 , 
                              bbp.550 = Fit.input$bbp.550, 
                              realdata = insitu.data, verbose = T,z = zB,
                              rb.fraction = fA.set)

rrs.forward.lee <- forward.op.lee[[1]]$Rrs #Extract Lee98 modelled Rrs

#4.5 Plot the forward simulated Rrs (AM03 & Lee98 with optional in situ)
if (insitu.present == TRUE) {
  forward.rrs <- data.frame("wave"=wavelength, "rrs.obs"=insitu.data,
                              "rrs.am03"=rrs.forward.am,
                              "rrs.lee98"=rrs.forward.lee)
  
  xmin = min(forward.rrs$wave); xmax= max(forward.rrs$wave); xstp=100
  ymin= 0; ymax=max(forward.rrs$rrs.am03)+0.20*max(forward.rrs$rrs.am03);ystp= signif(ymax/5, digits = 1)
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g1 <- ggplot()  + geom_line(data = forward.rrs,aes(y=rrs.obs,color="xx1",x = wave),
                              linetype="dashed",size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.am03,color="xx2",x = wave),
              size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.lee98,x = wave,color="xx3"), 
              size=1.3,show.legend = TRUE)+
    
    scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,actual"])),
                                   expression(paste(italic("R")["rs,model,AM03"])),
                                   expression(paste(italic("R")["rs,model,Lee98"]))), 
                        values = c("red","green","purple")) +
    
    scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+ 
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.position=c(0.55, 0.9),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_blank(),
          #legend.spacing.y = unit(2.0, 'cm'),
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  g1 
  
  if (plot == "TRUE") {
    ggsave(paste0("./Rrs_forward_chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
} else {
  forward.rrs <- data.frame("wave"=wavelength, #"rrs.obs"=insitu.data,
                            "rrs.am03"=rrs.forward.am,
                            "rrs.lee98"=rrs.forward.lee)
  
  xmin = min(forward.rrs$wave); xmax= max(forward.rrs$wave); xstp=100
  ymin= 0; ymax=max(forward.rrs$rrs.am03)+0.20*max(forward.rrs$rrs.am03);ystp= signif(ymax/5, digits = 1)
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g1 <- ggplot()  + 
    geom_line(data = forward.rrs,aes(y=rrs.am03,color="xx2",x = wave),
              size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.lee98,x = wave,color="xx3"), 
              size=1.3,show.legend = TRUE)+
    
    scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,AM03"])),
                                   expression(paste(italic("R")["rs,model,Lee98"]))), 
                        values = c("green","purple")) +
    
    scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+ 
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.position=c(0.55, 0.9),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_blank(),
          #legend.spacing.y = unit(2.0, 'cm'),
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  g1 
  
  if (plot == "TRUE") {
    ggsave(paste0("./Rrs_forward_chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
}

#------------------------------------------------------------------------------
#5. Perform inverse modeling using optimization of inverse cost function
#------------------------------------------------------------------------------

#5.1 Set Rrs observation data to invert
if (insitu.present == TRUE) {
  
  obsdata <-insitu.data #manual set observed rrs to begin inversion
  
}else {
  
  obsdata <-rrs.forward.am #set observed rrs to non-parametric SABER output
  
  #obsdata <-rrs.forward.am.sicf.dg_comp #set observed rrs to parametric SABER with SICF
}

#5.1 Pre-FIT of initial values
if (preFit == TRUE) {
  pre.Fit <- data.frame("C_ph"=seq(1,10,0.5), # <<USER DEFINED >>
                        "a_cdom.440"=seq(0.5,5,0.25),
                        "a.nap.440"=seq(0.01,0.1,0.005))
  
  pre.Fit.input.LUT  <- expand.grid(pre.Fit) #Create pre-Fit parameter space LUT
  
  preFIT.rrs.forward.LUT <- matrix(nrow = length(pre.Fit.input.LUT$C_ph),
                                   ncol = length(wavelength),0)
  #Create the Progress Bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(pre.Fit.input.LUT$C_ph), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")
  
  reslist = vector()
  for (i in 1:length(pre.Fit.input.LUT$C_ph)) { #Create Rrs LUT
    temp1 <- as.numeric(pre.Fit.input.LUT[i,])
    temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                           anap440 =temp1[3], bbp.550 = Fit.input$bbp.550,
                           realdata = obsdata,verbose=F )
    
    preFIT.rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
    reslist[i] = temp2[[2]]
    #cat(paste0("\033[0;43m",i," iterations over, ", (nrow(preFIT.rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
    setTxtProgressBar(pb, i)
    if (i == length(pre.Fit.input.LUT$C_ph)) {
      cat(paste0("\033[0;32m","###############PRE-FIT FINISHED################","\033[0m","\n"))
    }
  }
  
  prefit.best <- pre.Fit.input.LUT[which.min(reslist),] #retrieve best initial values 
                                                #using C.R.I.S.T.A.L.[minimizing SSR]
  prefit.best
  
  rrs.prefit <- Saber_forward(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440 , 
                              anap440 =prefit.best$a.nap.440, bbp.550 = Fit.input$bbp.550,
                              realdata = obsdata, verbose = T)[[1]]$Rrs
  
  #Show prefit spectra (Convert to ggplot2)
  plot(wavelength, obsdata, type="l", col="red", ylim=c(0,0.004))
  lines(wavelength, rrs.prefit, col="green")
}


#-----------------------------------------------------------------------------------------
##5.2 Prepare the optimization parameters
#-----------------------------------------------------------------------------------------

pop.sd = "unknown" ; constrain.inversion = TRUE 
manual_par0 = T

#5.2.1 Initial values for deep water 
if (pop.sd == "known" & type_Rrs_below == "deep" & manual_par0 == FALSE) { #@@ sigma KNOWN
  
  par0 = c(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440, 
           anap440 = prefit.best$a.nap.440)#, population.sd = 0.0006327431)
} 

if (pop.sd == "unknown" & type_Rrs_below == "deep" & manual_par0 == FALSE) { #@@ sigma UNKNOWN
  
  par0 = c(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440,
           anap440 = prefit.best$a.nap.440, population.sd = 0.001)
}

if (pop.sd == "unknown" & type_Rrs_below == "deep" & manual_par0 == TRUE) {# <<MANUAL INPUT>>
  par0 = c(chl = 4,
           acdom440 = 1, 
           anap440 =0.1,
           population.sd = 0.001) #@@  sigma UNKNOWN
}

#Set bounds (Instead of flat multiplier, implement range as a function 
                                            #of parameter sensitivity)
upper.bound <- par0 + 5*par0  
lower.bound <- par0 - 0.8*par0

#5.2.2 Initial values for shallow water
if (pop.sd == "unknown" & type_Rrs_below == "shallow") { # <<MANUAL INPUT>>
  par0 = c(chl = 2, acdom440 = 0.8,
           anap440 = 0.05, z = 2.5, 
           #rb.0 = 0.1,
           rb.1 = 0.5,
           rb.2 = 0.5,
           rb.3 = 0.5,
           rb.4 = 0.5,
           rb.5 = 0.5,
           population.sd = 0.1)
  
  #Autoscale Intital values from pre-Ft
  increament.scale <- 1
  
  lower.bound <- c((par0[1:3] - 0.8*par0[1:3]),z = 0.1,
                   #rb.0 = 0.1,
                   rb.1 = 0,
                   rb.2 = 0,
                   rb.3 = 0,
                   rb.4 = 0,
                   rb.5 = 0,
                   population.sd = 0.0001)
  
  upper.bound <- c((par0[1:3] + 5*par0[1:3]),z = 10,
                   #rb.0 = 1,
                   rb.1 = 1,
                   rb.2 = 1,
                   rb.3 = 1,
                   rb.4 = 1,
                   rb.5 = 1,
                   population.sd = 1)
}


#------------------------------------------------------------------------------------
#5.3 Gradient and finite difference based optimization 
#------------------------------------------------------------------------------------

obj = c("log-LL", "SSR"); obj.run <- obj[1]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt")

if (type_Rrs_below == "deep") {
  inverse_output <- suppressWarnings(solve.objective.inverse(initial = par0, 
                                            obsdata = obsdata,
                                            sa.model = "am03", obj.fn =obj.run , 
                                            method.opt = methods.opt[4],
                                            lower.b = lower.bound,
                                            upper.b = upper.bound, 
                                            batch = FALSE, pop.sd = FALSE))
}

if (type_Rrs_below == "shallow" & constrain.inversion == TRUE) {
  inverse_output <- solve.objective.inverse.shallow.constrained(constrained = T,
                                                                initial = as.numeric(par0), 
                                                    obsdata = obsdata,
                                                    sa.model = "lee98", obj.fn =obj.run , 
                                                    method.opt = methods.opt[4],
                                                    lower.b = lower.bound,
                                                    upper.b = upper.bound, 
                                                    batch = FALSE, pop.sd = FALSE)
}



if (obj.run == "log-LL" ) {
  if (type_Rrs_below == "deep") {
    Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$estimates[1], 
                                      "acdom.440"=inverse_output[[1]]$estimates[2],
                                      "anap.440"=inverse_output[[1]]$estimates[3])
  }
  if (type_Rrs_below == "shallow" & constrain.inversion == FALSE) {
    Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$estimates[1], 
                                      "acdom.440"=inverse_output[[1]]$estimates[2],
                                      "anap.440"=inverse_output[[1]]$estimates[3],
                                      "z"=inverse_output[[1]]$estimates[4],
                                      "rb.frac1"=inverse_output[[1]]$estimates[5],
                                      "rb.frac2"=inverse_output[[1]]$estimates[6],
                                      "rb.frac3"=inverse_output[[1]]$estimates[7],
                                      "rb.frac4"=inverse_output[[1]]$estimates[8],
                                      "rb.frac5"=inverse_output[[1]]$estimates[9]
                                      )
  }
  if (type_Rrs_below == "shallow" & constrain.inversion == TRUE) {
    Fit.optimized.ssobj <- data.frame("z"=inverse_output[[1]]$estimates[1],
                                      "rb.frac1"=inverse_output[[1]]$estimates[2],
                                      "rb.frac2"=inverse_output[[1]]$estimates[3],
                                      "rb.frac3"=inverse_output[[1]]$estimates[4],
                                      "rb.frac4"=inverse_output[[1]]$estimates[5],
                                      "rb.frac5"=inverse_output[[1]]$estimates[6]
    )
  
}} else{
  if (obj.run == "SSR") {
    
    Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$chl,
                                      "acdom.440"=inverse_output[[1]]$acdom.440,
                                "anap.440"=inverse_output[[1]]$anap.440)
  }
}

#------------------------------------------------------------------------------------
#5.4 Implement MCMC optimization
#------------------------------------------------------------------------------------
prior_fit = create.prior.data(use.ioccg.prior = F, use.wise.prior = F,
                  use.nomad.prior = T)

fit.chl.norm = prior_fit$fit.chl
fit.acdom440.norm = prior_fit$fit.acdom440
fit.anap440.norm = prior_fit$fit.anap440

#5.4.1 Create prior density and sampling class
if (pop.sd == "known" & type_Rrs_below == "deep") {
  prior.actual <- BayesianTools::createPrior(density = prior, sampler = sampler,
                                             lower = c(0,0,0), upper = c(30,5,0.5), 
                                             best = as.numeric(Fit.optimized.ssobj))
}

if (pop.sd == "unknown" & type_Rrs_below == "deep") {
  
  prior.actual <- BayesianTools::createPrior(density = prior, sampler = sampler,
                                             lower = c(0,0,0,0.0001),   # <<USER DEFINED>>
                                             upper = c(30,5,0.5, 0.01), # <<USER DEFINED>>
                                             #best = NULL)
                                             best = c(as.numeric(Fit.optimized.ssobj),
                                                      inverse_output[[1]]$estimates[4]))
}


#5.4.2 Create Bayesian setup for MCMC
#With prior
if (use.lklhood == FALSE) {
  bayessetup <- createBayesianSetup(prior = prior.actual,
                                    likelihood = ll,
                                    #lower = c(0,0,0), upper = c(30,5,0.5),
                                    names = c("chl","acdom440","anap440", "pop.sd"
                                              #, "x_not"
                                    ), 
                                    parallel = F)
}

#Only likelihood
if (use.lklhood == TRUE) {
  bayessetup <- createBayesianSetup(prior = NULL,
                                    likelihood = ll,
                                    lower = lower.bound, upper = upper.bound,
                                    names = c("chl","acdom440","anap440", "pop.sd"),
                                    parallel = F)
  
}

#5.4.3 Test if the setup is initiated for theta pars
checkBayesianSetup(bayessetup) 

#5.4.4 Set MCMC config
settings = list(iterations = 10000, message = TRUE, nrChains = 1, burnin=2000) 
samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

#5.4.5 Run MCMC
out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = samplerlist[6] )
summary(out)

#5.4.6 MCMC diagnostics
plot(out, start = 1000) #chain and parameter density
correlationPlot(out, start = 1000) #correlation plot among parameters
marginalPlot(out, start = 1000) #Variation in marginal prob density of prior and posterior

MAP.mcmc <- MAP(out) #Store MAP
DIC.mcmc <- DIC(out) #Store DIC

chain.one = as.data.frame(out[["codaChain"]][[1]])[-1,]
chain.two = as.data.frame(out[["codaChain"]][[2]])[-1,]
chain.three = as.data.frame(out[["codaChain"]][[3]])[-1,]

chain.mean = cbind(as.data.frame(out[["codaChain"]][[1]]),
                   as.data.frame(out[["codaChain"]][[2]]),
                   as.data.frame(out[["codaChain"]][[3]]))

chain = as.data.frame(t(apply(chain.mean,1, function(x) tapply(x,colnames(chain.mean),mean))))
chain = chain[-1,]

#------------------------------------------------------------------------------------------------
#Create MCMC plots (histograms)
#------------------------------------------------------------------------------------------------
png(filename = paste0("Posterior.SABER.chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), units = "in", width = 6.5, height = 4.5, res=300)
#dev.new()
par(mfrow = c(2,3))
#Posterior of [chl]
hist(chain$chl,nclass=30, main=expression(paste("Posterior of ", italic("[chl]")," [","mg"^1,"m"^-3,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$chl), col="purple", lwd=2)
abline(v = MAP(out)[[1]][1], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$chl, col="blue", lwd=2 )
abline(v = Fit.input$chl, col="black", lwd=2 )

#Posterior of acdom.440
hist(chain$acdom440,nclass=30, main=expression(paste("Posterior of ", italic(a)["CDOM"](440)," [","m"^-1,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$acdom440), col="purple", lwd=2)
abline(v = MAP(out)[[1]][2], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$acdom.440, col="blue", lwd=2 )
abline(v = Fit.input$acdom.440, col="black", lwd=2 )

#Posterior of anap.440
hist(chain$anap440,nclass=30, main=expression(paste("Posterior of ", italic(a)["NAP"](440)," [","m"^-1,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$anap440), col="purple", lwd=2)
abline(v = MAP(out)[[1]][3], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$anap.440, col="blue", lwd=2 )
abline(v = Fit.input$anap.440, col="black", lwd=2 )

#------------------------------------------------------------------------------------------------
#Create MCMC plots (chains)
#------------------------------------------------------------------------------------------------
#Chain of [chl]
plot(chain.one$chl, type = "l", xlab="iterations", col="orange", 
     #ylim=c(4.76,max(chain.mean$chl)),
     ylab= expression(paste(italic("[chl]"))),
     main = expression(paste("Chain of ", italic("[chl]")," [","mg"^1,"m"^-3,"]")))

lines(chain.two$chl, col="purple")
lines(chain.three$chl, col="grey")
abline(h = MAP(out)[[1]][1], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$chl, col="blue", lwd=2 )
abline(h = Fit.input$chl, col="black", lwd=2 )

#Chain of acdom.440
plot(chain.one$acdom440, type = "l", xlab="iterations", col="orange", ylab= expression(paste(italic(a)["CDOM"](440))),
     main = expression(paste("Chain of ", italic(a)["CDOM"](440)," [","m"^-1,"]")))
lines(chain.two$acdom440, col="purple")
lines(chain.three$acdom440, col="grey")
abline(h = MAP(out)[[1]][2], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$acdom.440, col="blue", lwd=2 )
abline(h = Fit.input$acdom.440, col="black", lwd=2 )

#Chain of anap.440
plot(chain.one$anap440, type = "l", xlab="iterations", col="orange", ylab= expression(paste(italic(a)["NAP"](440))),
     main = expression(paste("Chain of ", italic(a)["NAP"](440)," [","m"^-1,"]")))
lines(chain.two$anap440, col="purple")
lines(chain.three$anap440, col="grey")
abline(h = MAP(out)[[1]][3], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$anap.440, col="blue", lwd=2 )
abline(h = Fit.input$anap.440, col="black", lwd=2 )
dev.off()

fit.optimized.mcmc <- data.frame("chl"=MAP(out)[[1]][1], #Save MCMC MAP outputs
                                 "acdom.440"=MAP(out)[[1]][2],
                                 "anap.440"=MAP(out)[[1]][3])


#------------------------------------------------------------------------------------------
#6. Generate Rrs with the inversion retrieved parameters from both optimization
#------------------------------------------------------------------------------------------

#6.1 Generate Rrs with the inversion retrieved parameters 
if (type_Rrs_below == "deep") {
  #deep
  forward.fit.optimized.mle <- Saber_forward(chl = Fit.optimized.ssobj$chl, 
                                             acdom440 = Fit.optimized.ssobj$acdom.440, 
                                             anap440 =Fit.optimized.ssobj$anap.440,
                                             #bbp.550 = HL.deep.iop$bbp550[j],
                                             bbp.550 = Fit.input$bbp.550,
                                             realdata = obsdata)
}

if (type_Rrs_below == "shallow" & constrain.inversion == FALSE) {
  #Shallow
  forward.fit.optimized.mle <- Saber_forward(chl = Fit.optimized.ssobj$chl, 
                                             acdom440 = Fit.optimized.ssobj$acdom.440, 
                                             anap440 =Fit.optimized.ssobj$anap.440,
                                             #bbp.550 = HL.deep.iop$bbp550[j],
                                             z = Fit.optimized.ssobj$z,
                                             rb.fraction = as.numeric(Fit.optimized.ssobj[,5:9]),
                                             bbp.550 = Fit.input$bbp.550,
                                             realdata = obsdata, verbose = T, plot = T)
}

if (type_Rrs_below == "shallow" & constrain.inversion == TRUE) {
  #Shallow constrained
  forward.fit.optimized.mle <- Saber_forward(chl = Fit.input$chl, 
                                             acdom440 = Fit.input$acdom.440, 
                                             anap440 =Fit.input$anap.440,
                                             #bbp.550 = Fit.input$bbp.550,
                                             z = Fit.optimized.ssobj$z,
                                             rb.fraction = as.numeric(Fit.optimized.ssobj[-1]),
                                             bbp.550 = Fit.input$bbp.550,
                                             realdata = obsdata, verbose = T, plot = T)
}


rrs.forward.fit.optimized.mle <- forward.fit.optimized.mle[[1]]$Rrs


if (type_Rrs_below == "deep") {
  #deep
  forward.fit.optimized.mcmc <- Saber_forward(chl = fit.optimized.mcmc$chl, 
                                             acdom440 = fit.optimized.mcmc$acdom.440, 
                                             anap440 =fit.optimized.mcmc$anap.440,
                                             #bbp.550 = HL.deep.iop$bbp550[j],
                                             bbp.550 = Fit.input$bbp.550,
                                             realdata = obsdata)
}

if (type_Rrs_below == "shallow") {
  #Shallow
  forward.fit.optimized.mcmc <- Saber_forward(chl = fit.optimized.mcmc$chl, 
                                             acdom440 = fit.optimized.mcmc$acdom.440, 
                                             anap440 =fit.optimized.mcmc$anap.440,
                                             #bbp.550 = HL.deep.iop$bbp550[j],
                                             z = fit.optimized.mcmc$z,
                                             rb.fraction = as.numeric(fit.optimized.mcmc[,5:9]),
                                             bbp.550 = Fit.input$bbp.550,
                                             realdata = obsdata, verbose = T, plot = T)
}



rrs.forward.fit.optimized.mcmc <- forward.fit.optimized.mcmc[[1]]$Rrs


#6.2 Create Plot of actual Rrs and Rrs simulated with inversion retrieved parameters
#plot= FALSE
plotframe.rrs <- data.frame("wave"=wavelength, "rrs.est.sse"=rrs.forward.fit.optimized.mle,
                            "rrs.est.mcmc"=rrs.forward.fit.optimized.mcmc,
                            "rrs.obs"=obsdata, "rrs.prefit"= rrs.prefit)
#Create labels
mcmc.label = paste0("theta[MAP]== {",signif(Fit.optimized.mcmc$chl,digits = 2),"*',",signif(Fit.optimized.mcmc$acdom.440,digits = 2),",",signif(Fit.optimized.mcmc$anap.440,digits = 2),"'}")
mle.label = paste0("theta[MLE]== {",signif(Fit.optimized.ssobj$chl,digits = 2),"*',",signif(Fit.optimized.ssobj$acdom.440,digits = 2),",",signif(Fit.optimized.ssobj$anap.440,digits = 2),"'}")
obs.label = paste0("theta[obs]== {",signif(Fit.input$chl,digits = 2),"*',",signif(Fit.input$acdom.440,digits = 2),",",signif(Fit.input$anap.440,digits = 2),"'}")
prefit.label = paste0("theta[prefit]== {",signif(prefit.best$C_ph,digits = 2),"*',",signif(prefit.best$a_cdom.440,digits = 2),",",signif(prefit.best$a.nap.440,digits = 2),"'}")

#Create AXIS
xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
ymin= 0; ymax=max(plotframe.rrs$rrs.obs)+0.20*max(plotframe.rrs$rrs.obs)
ystp= signif(ymax/5, digits = 1)
asp_rat <- (xmax-xmin)/(ymax-ymin)

#Create Plot
g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est.sse,color="xx1",x = wave),
                            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.est.mcmc,color="xx2",x = wave),
            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.prefit,x = wave,color="xx5"), 
            size=1.3,show.legend = TRUE)+
  geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx6"),linetype="dashed", 
            size=1.3,show.legend = TRUE)+
  
  scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,MLE"])),
                                 expression(paste(italic("R")["rs,model,MAP"])),
                                 expression(paste(italic("R")["rs,model,prefit"])),
                                 expression(paste(italic("R")["rs,actual"]))), 
                      values = c("blue","green","purple","red")) +
  
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  
  annotate("text",x=470, y= ymax*0.95, label = obs.label, parse=T, size=4, color="red")+
  annotate("text",x=470, y= ymax*0.90, label = prefit.label, parse=T, size=4, color="purple")+
  annotate("text",x=470, y= ymax*0.85, label = mle.label, parse=T, size=4, color="blue")+
  annotate("text",x=470, y= ymax*0.80, label = mcmc.label, parse=T, size=4, color="green")+
  
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.62, 0.905),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        #legend.spacing.y = unit(2.0, 'cm'),
        plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1 

if (plot == "TRUE") {
  ggsave(paste0("./SABER.inversion.simulated.chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
}

#------------------------------------------------------------------------------------------
#7. Batch FORWARD RUN (TEST version, not for use)
#------------------------------------------------------------------------------------------
Fit <- data.frame("C_ph"=seq(1,10,0.5),
                  "a_cdom.440"=seq(0.5,5,0.25),
                  "a.nap.440"=seq(0.01,0.1,0.005))

Fit.input.LUT <- expand.grid(Fit) #Create Fit params LUT

rrs.forward.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(wavelength),0)

for (i in 1:length(Fit.input.LUT$C_ph)) { #Create Rrs LUT
  temp1 <- as.numeric(Fit.input.LUT[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3],bbp.550 = Fit.input$bbp.550, 
                         realdata = rrs.forward.am )
  rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
  
}

#------------------------------------------------------------------------------------------
##8. Batch Inversion (TEST version, not for use)
#------------------------------------------------------------------------------------------
Fit.optimized.ssobj.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(params),0)
colnames(Fit.optimized.ssobj.LUT) <- c("chl", "acdom440", "anap440")


for (i in 1:nrow(Fit.optimized.ssobj.LUT)) {
  
  inverse_output <- solve.objective.inverse(obj.fn = "log-LL", initial = par0, 
                                            obsdata =rrs.forward.LUT[i,] )
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(Fit.optimized.ssobj.LUT) - i), " remaining","\033[0m","\n"))
  Fit.optimized.ssobj.LUT[i,] <- as.numeric(inverse_output[[1]])
}


