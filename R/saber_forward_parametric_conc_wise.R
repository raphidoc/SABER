#=====================================================================================================
# Saber_forward_final.R simulates the remote sensing reflectance given the wavelengths,
# the water components along with bathymetry and bottom reflectance (for shallow water).

# This code has all possible modes "state-of-the-art" of bio-optical parametrizations for IOPs. It
# also has support to simulate Rrs using user provided in situ IOPs.

# The SA algorithm is extended for inelastic scattering equivalent Rrs as well. For SICF, Gilerson 2017
# is used whereas for fDOM, the code calculates analytical estimate. PLEASE NOTE, the analytical models
# are not depth integrated thus there can be possible discrepancies.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=====================================================================================================
#rm(list=ls(all=TRUE))
library(dplyr)
library(readxl)
library(stats4)
library(MASS)
library(dglm)
library(fitdistrplus)
library(Riops)
library(Cops)
library(ggplot2)
#--------------------------------------------------------------------------
#setwd("/home/musk0001/R_inverse_wasi")

# #Function to convert above water to under water geometry
# snell_law <- function(view,sun){
#   
#   # Index of refrations (real)
#   n_air= 1;  # air index of refration (real part)
#   n_w= 1.33; # water index of refration (real part)
#   
#   # Angles from the water
#   
#   # from deg to rad
#   view=view*(180/pi); # rad
#   sun=sun*(180/pi);   # rad
#   
#   # angles inside the water in rad
#   view_w= asin((n_air/n_w)*sin(view));   # rad
#   sun_w= asin((n_air/n_w)*sin(sun));     # rad
#   
#   # Fresnel Law
#   
#   rho_L= (1/2)*abs(((sin(view-view_w)^2)/(sin(view+view_w)^2))+((tan(view-view_w)^2)/(tan(view+view_w)^2)))
#   return(data.frame("view_w"=view_w, "sun_w"=sun_w, "rho_L"=rho_L))
# }

#Read default test data
demo.rrs = read.csv("./data/input-spectra/demo_rrs_Om.csv", header = T)

Saber_forward_final <-  function(use_true_IOPs = T, #Set TRUE if actual spectral IOPs exist
                                a_non_water_path = "./data/rb_retrieve_demo_a.csv", #a path
                                bb_non_water_path = "./data/rb_retrieve_demo_bb.csv", #bb path
                                           
                                           
                                chl=4.96, #must be input if SICF=TRUE and use_true_IOPs = F
                                dg_composite = TRUE, #Set TRUE if a_dg is provided together 
                                #when use_true_IOPs = F
                                           
                                acdom440=NULL, anap440=NULL, #if dg_composite = F, provide value of
                                # acdom440/443 (ag443) and anap440/443 (ad443)
                                
                                a_dg = 1,  #if dg_composite = T, provide value of
                                # acdom440/443 + anap440/443 (adg443)
                                
                                bbp.550=0.00726002, #bbp value at 550/555 nm
                                           
                                slope.parametric = TRUE, #Note, only possible if dg_composite = TRUE 
                                
                                use_spectral_shape_chl = FALSE, #Spectral shape normalization of a_phi
                                use_spectral_shape_dg = TRUE,   #Spectral shape normalization of a_dg
                                           
                                           
                                use_manual_slope = FALSE, #use manual spectral slopes
                                manual_slope = c("s_g"=0.015, "s_d"=0.01160, "gamma"=1), #Values of
                                #manual spectral slopes. must be provided in a named vector as shown
                                 
                                 
                                z=2, #bottom depth
                                rb.fraction = fA.set, #aerial fraction of bottom types
                                
                                use_spectral_rb = F, #use Spectral R_b instead of aerial fraction
                                spectral_rb_path = "./Outputs/Bottom_ref/Rb_spectral_data_MAN-R01.csv", #path
                                #for user-supplied bottom reflectance
                                
                                           
                                sicf = TRUE, q_phi=0.02, #sicf rrs and quantum yield
                                
                                fDOM = TRUE, #fdom rrs
                                
                                use_analytic_Ed = TRUE, #calculate Ed analytically from geometry
                                
                                sunzen_Ed = 60, #Sun Zenih at location, if avaiable, unless set -99
                                lat_Ed = 49, lon_Ed = -68, #Lat & Lon
                                date_time_Ed = "2019-08-18 20:50 GMT", #Date time in UTC
                        
                                Ed_fDOM_path = "./data/input-spectra/Ed_HL.csv", #Path for user supplied Ed file
                                use_fDOM_rad = F, #fDOM radiance is expected instead of reflectance
                                           
                                plot = FALSE, verbose = FALSE, #plot diagnostics and console output
                                
                                
                                realdata.exist = TRUE, #Set TRUE if in situ observation of rrs for which user wants
                                #to simulate, if such observation exists.
                                
                                realdata = demo.rrs$rrs.demo, # vector of realdata if realdata.exist = TRUE
                                # Note, The slope.parametric=T only works if valid data is 
                                # provided in realdata
                                
                                realdata_wave = seq(400,800,10) #wavelength of realdata if exists
                                
                                ){
  
  if (slope.parametric == TRUE & use_manual_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  
  ## OAC and IOP initialization
  
  # Read input params for forward modeling
  
  base.CDOM <- acdom440 # CDOM absorption at 440nm [m^-1]
  C_ph <- chl           # Phytoplankton concentration [mg/m^-3]
  base.NAP <- anap440   # NAP absorption at 440nm [m^-1]
  base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  lambda <- wavelength  #Desired spectral range for simulation
  
  if (verbose == T) {
    
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;36m","######################################### SIMULATION BEGINS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    
  }
    
  if (dg_composite == FALSE & slope.parametric == TRUE) {
    stop("QAA based paramteric derivation of absorption spectral slope can ONLY be calculated if dg_composite = TRUE")
  }
  
  
  if (use_true_IOPs == FALSE) {
    if (is.null(base.CDOM) & is.null(base.NAP)) {
      
      if (verbose == T) {
        print("Absorption of CDOM and NAP are read from a_dg argument")
      }
     
      base_CDM = a_dg
      
    } else {
      if (!is.null(a_dg)) {
        stop("Cannot run with values from both a_dg and acdom440, anap440; Set one set to NULL")
      } else {
        if (verbose == T) {
          print("absorption of CDOM and NAP are read from acdom440 and anap440 arguements")
        }
        
        base_CDM <- base.NAP + base.CDOM
      }
      
    }
  }
  
  ##If in situ IOPs are provided
  if (use_true_IOPs == TRUE) {
    
    if (!is.null(a_non_water_path) & !is.null(bb_non_water_path)) {
      
      a_non_water_df = read.csv(file = a_non_water_path, header = T)
      bb_non_water_df = read.csv(file = bb_non_water_path, header = T)
      
      if (is.data.frame(a_non_water_df) & is.data.frame(bb_non_water_df)) {
        
        a_non_water = a_non_water_df$at_w
        a_non_water_wave = a_non_water_df$wave
        
        bb_non_water = bb_non_water_df$bbp
        bb_non_water_wave = bb_non_water_df$wave
        
      } else {
        stop("The IOPs must be provided as data frames with <<wave>> as wavelength & <<a/bb>> as value")
      }
      
    }
  }
  
  #Shallow parameters
  fA <- rb.fraction     #Aerial fraction of bottom albedo
  zB <- z               #bottom depth
  
  #Actual AOP data for which we try to simulate the forward model
  if (!(all(length(realdata) == length(lambda)) && all(realdata_wave == lambda))) {
    
    if (verbose == TRUE) {
      print("Simulation wavelength and user-given wavelength of observed Rrs are different, data will be interpolated")
    }
    
    realdata_interp = Hmisc::approxExtrap(x= realdata_wave, y = realdata,
                                xout = lambda, method = "linear")$y
    
    realdata_interp[is.na(realdata_interp)] = 0
    Rrs_obs.interp <- realdata_interp 
    
  } else {
    
    Rrs_obs.interp <- realdata 
  }
  #lambda <- wavelength
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Pure water IOPs
  #--------------------------------------------------------------------------
  ## absorption (1/m)
  
  # wavelength range [190;4000] [nm]
  abs.water <- read.table("./data/input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  if (verbose == TRUE) {
    print(paste0("Water absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm.")) 
  }
  
  ## backscattering (1/m)
  
  if (type_case_water == 1) {# case1 water
    b1 <-  0.00144#  [1/m]
    
  } else {# case2 water
    if (type_case_water == 2) {
      b1 <-  0.00111#  [1/m]
    } else{
      print("Wrong water case type choosen")
    }
    
  }
  
  lambda1 <-  500# [nm]
  bb_W <- b1*(lambda/lambda1)^(-4.32)# [1/m]
  
  if (verbose == TRUE) {
    print(paste0("Water backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  #--------------------------------------------------------------------------
  ##Optically Active Constituents IOPs
  #--------------------------------------------------------------------------
  
  if (use_true_IOPs == FALSE) {
    
    if (verbose == TRUE) {
      print(paste0("Full spectral IOPs not provided, bio-optical models are used"))
    }
    
    ## Plankthon absorption (1/m)
    
    # load plankton absorption data
    A0_A1_PhytoPlanc <- read.table("./data/input-spectra/A0_A1_PhytoPlanc.dat")
    
    #Replace with the below for package environment for upper line
    #system.file("data", "input-spectra", "A0_A1_PhytoPlanc.dat", package = "SABER")
    
    # extract the values from the table
    lam_p <- A0_A1_PhytoPlanc$V1
    a0_p <- A0_A1_PhytoPlanc$V2
    a1_p <- A0_A1_PhytoPlanc$V3
    
    a0 <- rep(0,length(lambda))# [m^2/mg]
    a1 <- rep(0,length(lambda))# [m^2/mg]
    
    a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = lambda,method = "linear")$y
    a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = lambda,method = "linear")$y
    
    if (use_spectral_shape_chl == T) {
      
      aph_440 <- 1# [mg/m^3]
      abs_ph_norm <- rep(0,length(lambda))
      
      # Compute the value of plankton absorption as function of the concentration and wavelength
      for (j in 1:length(lambda)){
        
        abs_ph_norm[j] <- (a0[j] + a1[j]*log(aph_440))*aph_440
        
      }
      
      C_ph_norm=0.06*(C_ph)^0.65
      abs_ph = C_ph_norm*abs_ph_norm
      abs_ph[abs_ph < 0] <- 0
      
      if (verbose == T) {
        print("Planktonic absorption calculated using spectral shape scaled [chl]")
      }
      
      
    } else {
      
      # Compute the value of plankton absorption as function of the concentration and wavelength
      aph_440 <- 0.06*(C_ph)^0.65# [mg/m^3]
      abs_ph <- rep(0,length(lambda))
      
      for (i in 1:length(lambda)){
        
        abs_ph[i] <- (a0[i] + a1[i]*log(aph_440))*aph_440
        
      }
      abs_ph[abs_ph < 0] <- 0
      if (verbose == T) {
        print("Planktonic absorption calculated using absolute values of [chl]")
      }
      
    }
    
    
    if (verbose == TRUE) {
      print(paste0("Planktonic absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## CDOM+NAP absorption coefficient [1/m]
    
    if (dg_composite == TRUE) {
      if (use_spectral_shape_dg == TRUE) {
        
        Ga_CDOM <- base_CDM# [m^2/mg]
        Oa_CDOM <- 0
        
        if (slope.parametric == TRUE) {
          
          if (type_Rrs_below == "deep") {
            
            #parametric formula to retrieve spectral slope of CDOM + NAP
            S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
            
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
            }
            
            
          } else {
            if (verbose == T) {
              print("Shallow water type, QAA will derive wrong s_CDM, replaced with constant 0.017")
            }
            
            S_CDM <- 0.017 #<< HARD-CODED >>
          }
          
          
        } else {
          if (use_manual_slope  == TRUE) {
            S_CDM = as.numeric(manual_slope["s_g"] + manual_slope["s_d"])
            
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
            }
            
          } else {
            S_CDM <- 0.017 #<< Model Default >>
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
            }
            
          }
          
        }
        
        abs_CDM_440 <-  1 # [1/m], CDOM abs. coeff. at 440 [nm] normalized
        
        abs_CDM <- rep(0,length(lambda)); abs_CDM_norm <- rep(0,length(lambda))
        
        for (i in 1:length(lambda)){
          
          abs_CDM_norm[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
          
        }
        abs_CDM = Ga_CDOM*abs_CDM_norm
        
        if (verbose == T) {
          print(paste0("CDOM+NAP absorption calculated using spectral shape from ", min(lambda), " nm to ", max(lambda), " nm."))
        }
        
        
      } else {
        
        Ga_CDOM = 1
        Oa_CDOM = 0
        
        abs_CDM_440 <-  (Ga_CDOM*base_CDM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
        
        if (slope.parametric == TRUE) {
          
          if (type_Rrs_below == "deep") {
            #parametric formula to retrieve spectral slope of CDOM + NAP
            S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
            
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
            }
            
            
          } else {
            
            if (verbose == T) {
              print("Shallow water type, QAA will derive wrong s_CDM, replaced with constant 0.017")
            }
            
            S_CDM <- 0.017 #<< HARD CODED>>
            
          }
          
          
        } else {
          if (use_manual_slope  == TRUE) {
            S_CDM = as.numeric(manual_slope["s_g"] + manual_slope["s_d"])
            
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
            }
            
          } else {
            S_CDM <- 0.017 #<< Model Default >>
            
            if (verbose == T) {
              print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
            }
            
          }
          
        }
        
        abs_CDM <- rep(0,length(lambda))
        
        for (i in 1:length(lambda)){
          
          abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
          
        }
        
        if (verbose == T) {
          print(paste0("CDOM+NAP absorption calculated using absolute values from ", min(lambda), " nm to ", max(lambda), " nm."))
        }
        
        
      }
    } else {
      ## CDOM absorption coefficient [1/m]
      
      Ga_CDOM <- 1# [m^2/mg]
      Oa_CDOM <- 0
      
      if (use_manual_slope  == TRUE) {
        S_CDOM = as.numeric(manual_slope["s_g"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM is supplied by user as:", S_CDOM))
        }
        
      } else {
        S_CDOM <- 0.014 #<< Model Default >>
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM NAP is kept constant as:", S_CDOM))
        }
        
      }
      
      abs_CDOM_440 <-  (Ga_CDOM*base.CDOM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
      
      abs_CDOM <- rep(0,length(lambda))
      
      for (i in 1:length(lambda)){
        
        abs_CDOM[i]  <- abs_CDOM_440*exp(-S_CDOM*(lambda[i] - 440))
        
      }
      
      if (verbose == TRUE) {
        print(paste0("CDOM absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      ## SPM absorption coefficient [1/m]
      
      Ga_X <- 1# [m^2/mg]
      Oa_X <- 0
      
      
      if (use_manual_slope  == TRUE) {
        S_X = as.numeric(manual_slope["s_d"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for NAP is supplied by user as:", S_X))
        }
        
      } else {
        S_X <- 0.01160 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for NAP is kept constant as:", S_X))
        }
        
      }
      
      abs_X_440 <-  (Ga_X*base.NAP)+Oa_X# [1/m], SPM abs. coeff. at 440 [nm]
      
      abs_X <- rep(0,length(lambda))
      
      for (i in 1:length(lambda)){
        
        abs_X[i]  <- abs_X_440*exp(-S_X*(lambda[i] - 440))
        
      }
      
      if (verbose == TRUE) {
        print(paste0("NAP absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      print("Absorption of CDOM and NAP is calculated as individual components and constant slopes")
      abs_CDM = abs_CDOM + abs_X
    }
    
    abs_CDM_backup = abs_CDM #FOR fDOM
    
    #Plot the component specific absorption
    plotframe.abs <- data.frame("wave"=lambda,"abs.ph"= abs_ph, "abs.cdm"=abs_CDM)
    xmin = min(plotframe.abs$wave); xmax= max(plotframe.abs$wave); xstp=100
    ymin= 0; ymax=max(plotframe.abs$abs.cdm);ystp= signif(ymax/5, digits = 2)
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g <- ggplot()  + geom_line(data = plotframe.abs,aes(y=abs.ph,color="xx1",x = wave),
                               size=1.3,show.legend = TRUE) +
      geom_line(data = plotframe.abs,aes(y=abs.cdm,x = wave,color="xx2"), 
                size=1.3,show.legend = TRUE)+
      
      scale_colour_manual(labels = c(expression(paste(italic("a")[phi])),
                                     expression(paste(italic("a")["dg"]))
      ), 
      values = c("green","darkgoldenrod2")) +
      #ggtitle(paste0(stationlist[i])) +
      scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                         breaks = seq(xmin, xmax, xstp))  +
      scale_y_continuous(name =expression(paste(italic("a"),{}["t-w"],"(",lambda,",", 0^"-",")[", m^-1,"]")) , limits = c(ymin, ymax),
                         breaks = seq(ymin, ymax, ystp))+ 
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
            axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position=c(0.70, 0.9),
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
    g 
    
    if (plot == TRUE & dg_composite == FALSE) {
      ggsave(paste0("./Outputs/Forward/SABER.absorption.model_chl=",signif(C_ph, digits = 3) ,"_acdom440=", signif(base.CDOM, digits = 3),"_anap440=", signif(base.NAP,digits = 3), ".png"), plot = g,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
    }
    if (plot == TRUE & dg_composite == TRUE) {
      ggsave(paste0("./Outputs/Forward/SABER.absorption.model_chl=",signif(C_ph, digits = 3) ,"_adg440=", signif(base_CDM, digits = 3), ".png"), plot = g,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
    }
    
    #--------------------------------------------------------------------------
    ## Scattering and back-scattering
    #--------------------------------------------------------------------------
    
    # Backscattering coefficient for suspended particles [1/m]
    if (slope.parametric == TRUE) {
      
      if (type_Rrs_below == "deep") {
        #parametric formula to retrieve spectral slope of bbp using QAAv5
        refexponent = 2*(1-(1.2*exp(-0.9 * (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))]))))
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is calculated as:", refexponent))
        }
        
        
      } else {
        if (verbose == T) {
          print("Shallow water type, QAA will derive wrong spectral slope, manual slope will be used")
        }
        
        if (use_manual_slope  == TRUE) {
          refexponent = as.numeric(manual_slope["gamma"])
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
          }
          
        } else {
          refexponent <- 0.46 #<< Model Default >>
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
          }
          
        }
        
      }
      
    } else {
      
      if (use_manual_slope  == TRUE) {
        refexponent = as.numeric(manual_slope["gamma"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
        }
        
      } else {
        refexponent <- 0.46 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
        }
        
      }
    
    }
    
    
    bb_x <- rep(0,length(lambda))
    
    for (i in 1:length(lambda)){
      
      #print(lambda[i])
      bb_x[i] = base.bbp*((lambda[i]/550)^-(refexponent)) #Implement power model
      
    }
    if (verbose == TRUE) {
      print(paste0("NAP backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    #--------------------------------------------------------------------------
    ## Total Water system IOPs
    #--------------------------------------------------------------------------
    
    ## Total Absorption Coefficient (1/m)
    
    abs <-  a_W + abs_ph + abs_CDM #+ abs_X
    #abs= a_W + a_tw;
    
    if (verbose == TRUE) {
      print(paste0("Total absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## Total Backscattering Coefficient (1/m)
    
    bb <-  bb_W + bb_x
    #bb= bb_W + b_btw;
    
    if (verbose == TRUE) {
      print(paste0("Total backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
  } else {
    #--------------------------------------------------------------------------
    ## Water system total IOPs
    #--------------------------------------------------------------------------
    print(paste0("Actual IOPs provided"))
    
    a_non_water[a_non_water < 0] = 0
    
    ## Total Absorption Coefficient (1/m)
    if (!(all(length(a_non_water_wave) == length(lambda)) && all(a_non_water_wave == lambda))) {
      
      print("Simulation wavelength and absorption wavelength are different, data will be interpolated")
      
      a_non_water_interp = approx(x= a_non_water_wave, y = a_non_water,
                                  xout = lambda, method = "linear")$y
      
      a_non_water_interp[is.na(a_non_water_interp)] = 0
      abs <-  a_W + a_non_water_interp
      
    } else {
      
      abs <-  a_W + a_non_water
    }
    
    #Plot the fitted a_t-w
    if (plot == TRUE) {
      plot(a_non_water_df$wave, a_non_water_df$a, xlab="wavelength", 
           ylab="non-water absorption [m^-1]")
      lines(lambda, a_non_water_interp, col="red", lwd=3)
    }
    
    if (verbose == TRUE) {
      print(paste0("Total absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## Total Backscattering Coefficient (1/m)
    
    if ( !(all(length(bb_non_water_wave) == length(lambda)) && all(bb_non_water_wave == lambda))) {
      
      print("Simulation wavelength and backscatter wavelength are different")
      
      #### Compute bbp and bb spectral slope
      if (length(bb_non_water_wave) == 6) {
        print("HS-6 VSF is used")
        HS6_wl = c(394, 420, 470, 532, 620, 700)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      if (length(bb_non_water_wave) == 9) {
        
        print("BB9 VSF is used")
        HS6_wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      for (i in 1:nz) {
        #y = IOP.fitted.down$HS6$bbP[i,]
        y = as.matrix(bb_non_water)[,i]
        y[y < 0] = NA
        y[y > 0.1] = NA
        if (any(is.na(y))) {
          bbP555.down[i]=NA
          nuP.down[i]=NA
        } else {
          model = nls(y~b*x^z, start = list(b = y[3], z = 1),data=data.frame(x,y), 
                      control=list(maxiter=100, warnOnly=T))
          bbP555.down[i]=coef(model)[1]
          nuP.down[i]=coef(model)[2]
        }
      }
      #### Apply the power-law model to get hyperspectral bb
      waveletngth_hs <- matrix(data=lambda)
      refbbp <- matrix(data=bbP555.down)
      refexponent <- matrix(nuP.down)
      
      bbp_hs <- matrix(nrow = length(refbbp), ncol =length(waveletngth_hs),data = 0 )
      i=1
      j=1
      
      while (j<=length(refbbp)) {
        #print(refbbp[j,])
        for (i in 1:length(waveletngth_hs)) {
          #print(waveletngth_hs[i,])
          
          #Implement power model
          bbp_hs[j,i] = refbbp[j,1]*((waveletngth_hs[i,1]/555)^-(refexponent[j,1])) 
          
        }
        j=j+1
      }
      
      print("bbp for simulation wavelength is obtained using power-law fitting")
      bb <-  bb_W + as.vector(bbp_hs)
      
    } else {
      
      bb <-  bb_W + bb_non_water
    }
    
    #Plot the fitted bbp
    if (plot == TRUE) {
      plot(lambda, as.vector(bbp_hs),xlab="wavelength", 
           ylab="non-water backscatter [m^-1]")
      lines(bb_non_water_df$wave, bb_non_water_df$bb, col="red", lwd=3)
    }
    
    
    if (verbose == TRUE) {
      print(paste0("Total backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
  }
  
  ## Extinction Coefficient (1/m) and Single Back Scattering Albedo
  
  ext <-  abs+bb# [1/m] extinction coeff.
  omega_b <-  bb/ext# single back scattering albedo
  
  if (verbose == TRUE) {
    print(paste0("omega albedo calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ################################################################################################
  ##                                      RT Model                                              ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance below the surface
  #--------------------------------------------------------------------------
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w; view_w <- geometry$view_w; rho_L <- geometry$rho_L
  
  if (verbose == TRUE) {
    print("Viewing Geometry and Fresnel reflectance below surface calculated")
  }
  
  ## Remote Sensing Reflectance below the water surface
  
  if(type_case_water == "1"){
    #case 1
    
    f_rs <- 0.095 # [1/sr] 
  } else {
    if (type_case_water == "2") {
      #case 2
      
      f_rs <- 0.0512*(1 + (4.6659*omega_b) +(-7.8387*(omega_b^2)) + (5.4571*(omega_b^3)) )*(1+(0.1098/cos(sun_w)))*(1+(0.4021/cos(view_w)))# [1/sr]
    }
  }
  
  
  Rrs_below_deep <- f_rs*omega_b# [1/sr]
  
  if (type_Rrs_below == "deep") {
    
    Rrs_below <- Rrs_below_deep
    
    if (verbose == TRUE) {
      print(paste0("Subsurface (0^-)Rrs for deep water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
  } else{
    
    if (type_Rrs_below == "shallow") {
      
      #Bottom Contribution
      if (use_spectral_rb == FALSE) {
        
        if (verbose == T) {
          print("Bottom reflectance is assumed as LMM of user defined pure-spectra")
        }
        if (use_WASI_rb == TRUE) {
          if (verbose == T) {
            print("Bottom endmembers are taken from WASI default")
          }
          # Reflection factors of bottom surface [1/sr]
          #B0 <- 1/pi; 
          B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi; B4 <- 1/pi; B5 <- 1/pi; 
          #BOTTOM <- c(B0,B1,B2,B3,B4,B5)
          BOTTOM <- c(B1,B2,B3,B4,B5)
          
          # Bottom Albedo (costant)
          # wavelength range [350;900] [nm]
          bott0<- read.table("./data/input-spectra/Bott0const.R")
          wavebottom <- bott0$V1
          Bott0 <-  bott0$V2
          abott0 <- rep(0,length(lambda))
          abott0 <-  Hmisc::approxExtrap(wavebottom, Bott0, xout = lambda, method = "linear")$y
          
          # Bottom Albedo Sand
          # wavelenght range [350;1000] [nm]
          bott1 <- read.table("./data/input-spectra/Bott1SAND.R")
          wavebottom <- bott1$V1
          Bott1 <-  bott1$V2
          abott1 <- rep(0,length(lambda))
          abott1 <-  Hmisc::approxExtrap(wavebottom, Bott1, xout = lambda, method = "linear")$y
          
          # Bottom Albedo of fine-grained sediment
          # wavelenght range [350;900] [nm]
          bott2 <- read.table("./data/input-spectra/Bott2silt.R")
          wavebottom <- bott2$V1
          Bott2 <-  bott2$V2
          abott2 <- rep(0,length(lambda))
          abott2 <-  Hmisc::approxExtrap(wavebottom, Bott2, xout = lambda, method = "linear")$y
          
          # Bottom Albedo of green makrophyte "Chara contraria"
          # wavelenght range [350;900] [nm]
          bott3 <- read.table("./data/input-spectra/Bott3chara.R")
          wavebottom <- bott3$V1
          Bott3 <-  bott3$V2
          abott3 <- rep(0,length(lambda))
          abott3 <-  Hmisc::approxExtrap(wavebottom, Bott3, xout = lambda, method = "linear")$y
          
          # Bottom Albedo of green makrophyte "Potamogeton perfoliatus"
          # wavelenght range [350;900] [nm]
          bott4 <- read.table("./data/input-spectra/Bott4perfol.R")
          wavebottom <- bott4$V1
          Bott4 <-  bott4$V2
          abott4 <- rep(0,length(lambda))
          abott4 <-  Hmisc::approxExtrap(wavebottom, Bott4, xout = lambda, method = "linear")$y
          
          # Bottom Albedo of green makrophyte "Potamogeton pectinatus"
          # wavelenght range [350;900] [nm]
          bott5 <- read.table("./data/input-spectra/Bott5pectin.R")
          wavebottom <- bott5$V1
          Bott5 <-  bott5$V2
          abott5 <- rep(0,length(lambda))
          abott5 <-  Hmisc::approxExtrap(wavebottom, Bott5, xout = lambda, method = "linear")$y
          
          #abott <- rbind(abott0, abott1, abott2, abott3, abott4, abott5)
          abott <- rbind(abott1, abott2, abott3, abott4, abott5)
          
          Bottom <-  matrix(nrow = nrow(abott), ncol = ncol(abott), 0)
          Rrs_Bottom <- matrix(nrow = nrow(abott), ncol = ncol(abott), 0)# Bottom remote sensing reflectance [1/sr]
          
          for (i in 1:length(fA)){
            Bottom[i,] <-  fA[i]*abott[i,]
            Rrs_Bottom[i,] <-  BOTTOM[i]* Bottom[i,] #fA(i)*abott(:,i);
          }
          
          Bottom <- colSums(Bottom)
          Rrs_Bottom <- colSums(Rrs_Bottom)# [1/sr]
          
        } else {
          if (verbose == T) {
            print("Bottom endmembers are taken from WISE-Man/Algae-WISE")
          }
          
          # Reflection factors of bottom surface [1/sr]
          #B0 <- 1/pi; 
          B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi #;B4 <- 1/pi; B5 <- 1/pi; 
          BOTTOM <- c(B1,B2,B3)#,B4,B5)
          
          # Bottom Albedo for ALGAE-WISE
          
          abott1 <-  rb$class1 
          abott2 <-  rb$class2
          abott3 <-  rb$class3 
          
          abott <- rbind(abott1, abott2, abott3)#, abott4, abott5)
          
          Bottom <-  matrix(nrow = length(fA), ncol = ncol(abott), 0)
          Rrs_Bottom <- matrix(nrow = length(fA), ncol = ncol(abott), 0)# Bottom Rrs [1/sr]
          
          Bottom = fA * abott #BitWISE operation
          
          Rrs_Bottom = BOTTOM * Bottom #BitWISE operation
          
          Bottom <- colSums(Bottom) #[unitless]
          Rrs_Bottom <- colSums(Rrs_Bottom)# [1/sr]
          
        }
    
        
        
      } else {
        
        if (verbose == T) {
          print("Bottom reflectance is read from user-defined bottom reflectance file")
        }
        
        rb_bottom = read.csv(spectral_rb_path, header = T)
        rb_bottom_wave = rb_bottom$wave
        rb_bottom_rb = rb_bottom$rb_est_iop
        
        if ( !(all(length(rb_bottom_wave) == length(lambda)) && all(rb_bottom_wave == lambda))) {
          print("Simulation wavelength and supplied bottom reflectancewavelength are different")
          rb_bottom_rb_interp = approx(x=rb_bottom_wave, y = rb_bottom_rb, xout = lambda, 
                                       method = "linear")$y
          Rrs_Bottom = rb_bottom_rb_interp
        } else {
          Rrs_Bottom = rb_bottom_rb
        }
        
      }
      
      if (verbose == TRUE) {
        print(paste0("Aerial fraction scaled bottom albedo generated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      # Attenuation Coefficients
      if (type_case_water == 1) {
        
        k0 <- 1.0395 #case 1
      } else {
        if (type_case_water == 2) {
          
          k0 <- 1.0546 #case 2
        }
      }
      
      Kd <- k0*(ext/cos(sun_w))
      kuW <- (ext/cos(view_w))*((1+omega_b)^3.5421)*(1-(0.2786/cos(sun_w)))
      kuB <- (ext/cos(view_w))*((1+omega_b)^2.2658)*(1-(0.0577/cos(sun_w)))
      
      if (verbose == TRUE) {
        print(paste0("Attenuation coefficients calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      
      #Final calculation for shallow Rrs
      Ars1 <- 1.1576; Ars2 <- 1.0389 #Paramteric coeffs for shallow water
      Rrs_below_shallow <-  Rrs_below_deep*(1-(Ars1*exp(-zB*(Kd+kuW)))) + Ars2*Rrs_Bottom*exp(-zB*(Kd+kuB)) 
      Rrs_below <- Rrs_below_shallow
      
      if (verbose == TRUE) {
        print(paste0("Subsurface (0^-) Rrs in Shallow water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
    }
  }
  
  #Store the Rrs from elastic scattering
  Rrs_elastic = Rrs_below
  
  #===============================================================
  #Calculate [chl] Fluorescence equivalent Rrs
  #===============================================================
  if (sicf == TRUE & dg_composite == FALSE) {
    if (is.null(C_ph)) {
      stop("Chl concentration must be provided to calculate SICF")
    }
    rrs_sicf = Rrs_Fluorescence(c_chl = C_ph,  
                                #c_chl = Fit.input$chl,
                                dg_comsposite = FALSE,  wavelength = lambda,
                                dg_443 = NULL,
                                
                                use_analytic_Ed = TRUE,
                                sunzen_Ed = sunzen_Ed, lat_Ed = lat_Ed, lon_Ed = lon_Ed,
                                date_time_Ed = date_time_Ed,
                                
                                abs_cdom_443 = base.CDOM,
                                abs_nap_443 = base.NAP, phi_f = q_phi)
    if (verbose == T) {
      print(paste0("SICF: Absorption of CDOM and NAP is user defined inputs"))
      print(paste0("Subsurface (0^-) Rrs equivalent to SICF calculated with quantum yield of ", q_phi))
    }
    
    Rrs_below = Rrs_below + rrs_sicf
  }
  
  if (sicf == TRUE & dg_composite == TRUE) {
    if (is.null(C_ph)) {
      stop("Chl concentration must be provided to calculate SICF")
    }
    qaa_op = QAA.v5(waves = lambda, Rrs = Rrs_obs.interp)
    rrs_sicf = Rrs_Fluorescence(c_chl = C_ph,  
                                #c_chl = Fit.input$chl,
                                dg_comsposite = TRUE,  wavelength = lambda,
                                
                                use_analytic_Ed = TRUE,
                                sunzen_Ed = sunzen_Ed, lat_Ed = lat_Ed, lon_Ed = lon_Ed,
                                date_time_Ed = date_time_Ed,
                                
                                dg_443 = qaa_op$a_dg_443,
                                abs_cdom_443 = base.CDOM,
                                abs_nap_443 = base.NAP, phi_f = q_phi)
    if (verbose == T) {
      print(paste0("SICF: Absorption of CDOM and NAP is calculated from QAA"))
      print(paste0("Subsurface (0^-) Rrs equivalent to SICF calculated with quantum yield of ", q_phi))
    }
    
    Rrs_below = Rrs_below + rrs_sicf
  }
  
  #===============================================================
  #Calculate CDOM Fluorescence equivalent Rrs
  #===============================================================
  if (fDOM == TRUE) {
    print("Calculating the discretized wavelength redistribution function...")
    #Implement the wavelength redistribution functions with discretization 
    #for inelastic scattering 
    wave_redist_inel = wave.redist.inelastic(quant_y_phi = q_phi, include.chl = T,
                                             include.CDOM = T, wavelength = lambda)
    
    #Plot the discretized wavelength redistribution functions for inelastic scattering 
    if (plot == TRUE) {
      plot(lambda[-length(lambda)], wave_redist_inel[["sicf"]][1,], pch=19, type="l", 
           lwd=2.5, col="green",
           ylab = expression(paste("Fluoroscence Efficieny (1/nm)")), xlab = "Wavelength" )
      
      lines(lambda[-length(lambda)], wave_redist_inel[["fDOM"]][1,], col="goldenrod2", lwd=2.5)
      
      legend(x = "topleft",          # Position
             legend = c("SICF","fDOM"),  # Legend texts
             lty = c(1, 1),           # Line types
             col = rev(c("goldenrod2", "green")),           # Line colors
             lwd = 2)
    }
    
    #=====================================================================
    #Obtain the fDOM scattering function components
    #=====================================================================

    waveb <-lambda
    Nwave <- length(waveb) - 1
    
    WRF_CDOM = wave_redist_inel$fDOM # extract the wave discretized quantum yield shape for fDOM
    
    ###### Scalar irradiance calculation #####
    
    #---------------------------------------------------------------------------
    #Get Ed0 from HL simulated Ed data
    #---------------------------------------------------------------------------
    if (use_analytic_Ed == FALSE) {
      
      Ed_path = Ed_fDOM_path
      
      print("NOTE: Pre-Defined Ed is used; results might be erroneous")
      
      Ed_sim = read.csv(file = Ed_path, header = TRUE, skip = 9)
      Ed_sim_interp = Hmisc::approxExtrap(x = Ed_sim$Wavelength, y = log(Ed_sim$Ed_total.W.m.2.nm.),
                                          xout = wavelength, method = "linear")$y
      Ed_sim_interp = exp(Ed_sim_interp)*100 #extrapolate and change units
      
      Ed_sim_interp_0m = 0.96*Ed_sim_interp #Convert to underwater E0
      
      sun_view = snell_law(sun = sunzen_below(sun_zen_aove = sun_above), view = view) #calculate sun zenith
      
      E0_sim_interp_scalar = Ed_sim_interp_0m/(cos(sun_view$sun_w)) #Convert to scalar irradiance
      
      E0_0m = E0_sim_interp_scalar
      Ed_0m = Ed_sim_interp_0m
      
    } else {
      
      #---------------------------------------------------------------------------
      #Calculate the Ed using Gergg & Carder model
      #---------------------------------------------------------------------------
      
      print("NOTE: Gregg & Carder 1990 is used to estimate Ed")
      
      Cops::GreggCarder.data()
      
      library(lubridate)
      library(dplyr)
      
      x = as_datetime(as.POSIXct(date_time_Ed))
      
      jday_no = yday(x) #Get Jullien Day no.
      time_no = format(x, "%T") #Get time in UTC
      
      time_dec = sapply(strsplit(as.character(time_no),":"), #Convert time in decimal format
             function(x) {
               x <- as.numeric(x)
               x[1]+x[2]/60
             }
      )
      
      test_Ed = GreggCarder.f.modified(the = sunzen_Ed, #Calculate the Ed following Gregg & Carder 1990
                              lam.sel = lambda, hr = time_dec,
                              jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)

      if (sunzen_Ed < 0) {
        print("Sun Zenith was not provided, calculated from Geometry")
        sunzen_Ed = Cops::GreggCarder.sunang(rad = 180/pi, iday = jday_no, 
                                             xlon = lon_Ed, ylat = lat_Ed, hr = time_dec)
      }
      Ed0.0p = test_Ed$Ed #Total Ed at 0+
      Ed0_dir.0p = test_Ed$Edir #Direct Ed at 0+
      Ed0_dif.0p = test_Ed$Edif #Diffused Ed at 0+
      
      rhoF <- GreggCarder.sfcrfl(rad = 180.0/pi, theta=sunzen_Ed, ws=5) #Calculate fresnel reflectance
      sun_view_op = snell_law(view = 0 , sun = sunzen_Ed) #calculate sun angle 0+ in radian
      
      #Translate Ed0+ to Ed0-
      Ed0.0m = (Ed0_dir.0p *(1 - rhoF$rod)) + # direct
        (Ed0_dif.0p * (1 - rhoF$ros))   # diffuse
      
      #Convert it to scalar irradiance
      sun_view_om = snell_law(view = 0 , sun = sunzen_below(sunzen_Ed))
      Ed0m_0 = (Ed0.0m*100) / cos(sun_view_om$sun_w)
      
      #Plot the simulated Ed profiles
      if (plot == TRUE) {
        plot(lambda, Ed0.0p*100, col="black", pch=19, ylab = "Ed", ylim=c(0,1.25*max(Ed0m_0)))
        lines(lambda, Ed0.0m*100, col="red", lwd=2.5, lty="dashed")
        lines(lambda, Ed0m_0, col="navyblue", lwd=2.5)
        
        legend(x = "bottom",          # Position
               legend = c("Ed0+","Ed0-", "E0"),  # Legend texts
               lty = c(1, 2, 1),           # Line types
               col = c("black", "red", "navyblue"), # Line colors
               lwd = 2, bty = "n")
        
      }
      
      E0_0m = Ed0m_0
      Ed_0m = Ed0.0m*100

    }
    
    ###### absorption of CDOM calculation #####
    
    if ( use_true_IOPs == TRUE ) {
      qaa_op = QAA.v5(waves = lambda, Rrs = Rrs_obs.interp)
      abs_CDM_440 =  qaa_op$a_dg_443
      
      if (type_Rrs_below == "deep") {
        if (slope.parametric == TRUE) {
          #parametric formula to retrieve spectral slope of CDOM + NAP
          S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
          print(paste0("fDOM: The spectral slope for CDOM + NAP is calculated as:", S_CDM))
        } else{
          
          if (use_manual_slope  == TRUE) {
            S_CDM = as.numeric(manual_slope["s_g"] + manual_slope["s_d"])
            print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
          } else {
            S_CDM <- 0.017 #<< Model Default >>
            print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
          }
          
        }
        
        
      } else {
        print("fDOM: The water type is shallow, QAA will obtain wrong s_CDOM, thus const. 0.017 is used")
        S_CDM <- 0.017 #<< USER INPUT >>
      }
      
      abs_CDM <- rep(0,length(lambda))
      
      for (i in 1:length(lambda)){
        
        abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
        
      }
      
    } else {
      #abs_CDM <- rep(0,length(lambda))
      abs_CDM = abs_CDM_backup
      
    }
    
    sum0 = matrix(0, ncol = Nwave, nrow = Nwave) #stores the lambda specific integration of fDOM 
    
    #######Do the integration with respect to wavelength of the final source function of fDOM#####
    # the equation to be integrated is:a_CDOM*f_F(ex, em)
    
    for (jwave in 2:Nwave) { #loop for excitation (lambda_ex)
      
      for (iwave in 1:(jwave-1)) { #loop for emission (lambda_exm)
        
        Eozi <- E0_0m[iwave] # retrieve E0 at lambda_ex
        
        absCDOM <- abs_CDM[iwave] # retrieve a_CDOM at lambda_ex
        
        sum0[iwave,jwave] <- sum0[iwave,jwave] + (Eozi * absCDOM * WRF_CDOM[iwave, jwave]) #Integrate
        
      }
    }
    
    fDOM_rad = sum0[1,]/(4*pi) #multiply with fluoroscence anisotropic phase function
    print("fDOM equivalent radiance is calculated")
    
    #fDOM_rrs = fDOM_rad / Ed_interp[-length(Ed_interp)]
    #fDOM_rrs = fDOM_rad / Ed0_m_interp[-length(Ed0_m_interp)]
    fDOM_rrs = fDOM_rad / Ed_0m[-length(Ed_0m)]
    
    if (use_fDOM_rad == FALSE) {
      
      #Rrs_below = Rrs_below + fDOM_rrs
      
      rrs_final = Rrs_below[1:length(Rrs_below)-1] + fDOM_rrs
      
      # output_rrs = fDOM_rrs
      # output_rrs[-length(output_rrs)] = rrs_final
      
      Rrs_below[1:length(Rrs_below)-1] = rrs_final
      
    } else {
      
      #Rrs_below = Rrs_below + fDOM_rad
      
      rrs_final = Rrs_below[1:length(Rrs_below)-1] + fDOM_rad
      
      # output_rrs = fDOM_rad
      # output_rrs[-length(output_rrs)] = rrs_final
      
      Rrs_below[1:length(Rrs_below)-1] = rrs_final
      
      
    }
    
    print(paste0("Subsurface (0^-) Rrs equivalent to fDOM is calculated"))
    
    
  }
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance above the surface
  #--------------------------------------------------------------------------
  # Extraterrestrial solar irradiance [mW/m^2 nm]
  E0 <-  read.table("./data/input-spectra/E0.txt", header = F)
  E0.wavelength <-  E0$V1
  E0.ett <-  E0$V2
  E0 <- rep(length(lambda), 0)
  
  E0 <-  Hmisc::approxExtrap(x =E0.wavelength, y =E0.ett, xout = lambda, method = "linear")$y
  
  # Oxygen absorption [1/cm]
  absO2 <- read.table("./data/input-spectra/absO2.A", header = F)
  absO2.wavelength <-  absO2$V1
  absO2.oxy <-  absO2$V2
  abs_O2 <- rep(length(lambda), 0)
  
  abs_O2 <-  Hmisc::approxExtrap(x =absO2.wavelength, y =absO2.oxy, xout = lambda, method = "linear")$y
  
  # Ozone absorption [1/cm]
  absO3 <- read.table("./data/input-spectra/absO3.A", header = F)
  absO3.wavelength <-  absO3$V1
  absO3.oxy <-  absO3$V2
  abs_O3 <- rep(length(lambda), 0)
  
  abs_O3 <-  Hmisc::approxExtrap(x =absO3.wavelength, y =absO3.oxy, xout = lambda, method = "linear")$y
  
  # Water vapour absorption [1/cm]
  absWV <- read.table("./data/input-spectra/absWV.A", header = F)
  absWV.wavelength <-  absWV$V1
  absWV.wv <-  absWV$V2
  abs_WV <- rep(length(lambda), 0)
  
  abs_WV <-  Hmisc::approxExtrap(x =absWV.wavelength, y =absWV.wv, xout = lambda, method = "linear")$y
  
  # angles from deg to rad
  
  view_rad=view*(180/pi); # rad
  sun_rad=sun*(180/pi);   # rad
  
  # Downwelling Irradiance [mW/m^2 nm]
  M <- 1/(cos(sun_rad)+(0.50572*((90+ 6.079975-sun)^(-1.253))))
  M1 <-  (M*P)/1013.25
  Moz <-  1.0035/ (((cos(sun_rad)^2)+0.007)^0.5)
  
  
  # Air mass type
  if (type_case_water == "1") { 
    AM <- 1
  } else{
    
    AM <- 1
  }
  
  omega_a <- ((-0.0032*AM) + 0.972)*exp(RH*3.06*(10^-4))
  
  Ha <- 1; V <- 15# [km] aerosol scale height and horizontal visibility
  beta <- 3.91*(Ha/V)
  tau_a <-  beta*((lambda/550)^(-alpha))
  
  Tr <-   exp(-M1/((115.6406*(lambda^(4)))-(1.335*(lambda^(2)))))
  Taa <-  exp(-(1-omega_a)*tau_a*M)
  Tas <-  exp(-omega_a*tau_a*M)
  Toz <-  exp(-abs_O3*Hoz*Moz)
  To <-   exp((-1.41*abs_O2*M1)/((1+(118.3*abs_O2*M1))^0.45))
  Twv <-  exp((-0.2385*abs_WV*WV*M)/((1+(20.07*abs_WV*WV*M))^0.45))
  
  
  B3 <- 0.82 - (0.1417*alpha); B1 <-  B3*(1.459 +(B3*(0.1595+(0.4129*B3)))); B2 <-  B3*(0.0783 +(B3*(-0.3824-(0.5874*B3))))
  Fa <-  1-(0.5*exp((B1+(B2*cos(sun_rad)))*cos(sun_rad)))
  
  Edd <-  E0*Tr*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsr <- (1/2)*E0*(1-(Tr^(0.95)))*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsa <-  E0*(Tr^(1/2))*Taa*(1-Tas)*Toz*To*Twv*cos(sun_rad)*Fa
  Eds <-  Edsr+Edsa # diffuse downwelling irradiance (sum of Rayleigh and aerosol)
  
  Ed <-  (f_dd*Edd) + (f_ds*Eds)# downwelling irradiance
  
  # Sky Radiance
  Ls <-  (g_dd*Edd) + (g_dsr*Edsr) + (g_dsa*Edsa)# [mW/sr m^2 nm]
  
  # Remote sensing reflectance above the surface [1/sr]
  Rrs_above <-  rho_L*(Ls/Ed)
  
  if (verbose == TRUE) {
    print(paste0("Skyglint calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------
  sigma <- 0.03; nW <- 1.33; rho_U <- 0.54; Q <- 5#[sr]
  
  if (type_Rrs_water == "below_surface") {
    
    Rrs <-  Rrs_below
    
    if (verbose == TRUE) {
      print(paste0("FINAL: Subsurface (0^-)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
  } else {
    if (type_Rrs_water == "above_surface_with_glint") {
      
      Rrs <-  (((1-sigma)*(1-rho_L)/(nW^2))* (Rrs_below/(1-rho_U*Q*Rrs_below)))+ Rrs_above
      
      if (verbose == TRUE) {
        print(paste0("FINAL: Above surface (0^+)Rrs with modelled glint calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
    } else {
      if (type_Rrs_water == "above_surface_only") {
        
        Rrs <-  (((1-sigma)*(1-rho_L)/(nW^2))* (Rrs_below/(1-rho_U*Q*Rrs_below)))
        
        if (verbose == TRUE) {
          print(paste0("FINAL: Above surface (0^+)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
          
        }
      }
    }
  }
  if (verbose == T) {
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;32m","######################################### SIMULATION ENDS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
  }
  
  
  #Compare actual vs modelled Rrs
  #Rrs_obs.interp <- Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs, xout = lambda, method = "linear")$y
  if (realdata.exist == TRUE) {
    #browser()
    plotframe.rrs <- data.frame("wave"=lambda, "rrs.est"=Rrs, "rrs.obs"=Rrs_obs.interp)
    xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
    ymin= 0; ymax=max(plotframe.rrs$rrs.obs) 
    ymax=signif(ymax + 0.5*ymax, digits = 1)
    ystp= signif(ymax/5, digits = 1)
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est,color="xx1",x = wave),
                                size=1.3,show.legend = TRUE) +
      geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx5"),linetype="dashed", 
                size=1.3,show.legend = TRUE)+
      scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model"])),
                                     expression(paste(italic("R")["rs,actual"]))), 
                          values = c("blue","green")) +
      #ggtitle(paste0(stationlist[i])) +
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
            legend.position=c(0.70, 0.9),
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
    if (plot == TRUE & use_true_IOPs == TRUE) {
      ggsave(paste0("./Outputs/Forward/SABER.forward.Rrs_station_",statname, ".png"), plot = g1,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
      
    } else {
      if (plot == "TRUE" & type_Rrs_below == "deep") {
        ggsave(paste0("./Outputs/Forward/SABER.forward.Rrs_chl_",signif(C_ph, digits = 3) ,"_adg440_", signif(base_CDM, digits = 3),".png"), plot = g1,
               scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
      }
      if (plot == "TRUE" & type_Rrs_below == "shallow") {
        ggsave(paste0("./Outputs/Forward/SABER.forward.Rrs_chl_",signif(C_ph, digits = 3) ,"_adg440_", signif(base_CDM, digits = 3), "_z_", signif(zB,digits = 3),  ".png"), plot = g1,
               scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
      }
      
    }
    
    #--------------------------------------------------------------------------
    ## Objective function for the inverse mode
    #--------------------------------------------------------------------------
    #Calculate minimizable residual function
    Res <- sum((Rrs_obs.interp-Rrs)^2) #euclidean space
    
    #Calculate residual as function of wavelength
    Res.spectral <- (Rrs - Rrs_obs.interp)*100/Rrs_obs.interp
    
    if (use_true_IOPs == F) {
      return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs, "Rrs_elastic" = Rrs_elastic,
                             "p.bias"=Res.spectral),"abs_comp"=plotframe.abs, 
                  "bbp"= bb_x, 
                  "ss.residual"=Res,"method"=c("SSR eucledian")))
    } else {
      return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs, "Rrs_elastic" = Rrs_elastic,
                             "p.bias"=Res.spectral),
                  "a_t" = (abs - a_W), "bbp" = (bb - bb_W),
                  "ss.residual"=Res,"method"=c("SSR eucledian")))
    }
    
    #return(Res)
    
  } else {
    plotframe.rrs <- data.frame("wave"=lambda, "rrs.est"=Rrs)
    xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
    ymin= 0; ymax=max(plotframe.rrs$rrs.est);ystp= signif(ymax/5, digits = 1)
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est,color="xx1",x = wave),
                                size=1.3,show.legend = TRUE) +
      # geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx5"),linetype="dashed", 
      #           size=1.3,show.legend = TRUE)+
      scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model"]))),
                          #expression(paste(italic("R")["rs,actual"]))), 
                          values = c("blue"))+
      #,"green")) +
      #ggtitle(paste0(stationlist[i])) +
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
            legend.position=c(0.70, 0.9),
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
    
    if (plot == "TRUE" & type_Rrs_below == "deep") {
      ggsave(paste0("./Outputs/Forward/SABER.forward.Rrs_chl_",signif(C_ph, digits = 3) ,"_adg440_", signif(base_CDM, digits = 3), ".png"), plot = g1,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
    }
    if (plot == "TRUE" & type_Rrs_below == "shallow") {
      ggsave(paste0("./Outputs/Forward/SABER.forward.Rrs_chl_",signif(C_ph, digits = 3) ,"_adg440_", signif(base_CDM, digits = 3), "_z_", signif(zB,digits = 3),  ".png"), plot = g1,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
    }
    
    return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs)))
    
  }
  
  
}
