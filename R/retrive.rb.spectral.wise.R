#==================================================================================================
# retrieve.rb.spectral.wise.R performs spectral retrieval of bottom reflectance (Rb) based on 
# Radiative transfer in optically complex case-II waters. There are two modes to retrieve Rb:

# 1. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given [chl],
# CDOM and NAP absorption at 440nm wavelength and depth and subtract the model value from observed 
# Rrs using AM03. The spectral slopes can either be parametric or user provided constants. 

# 2. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given spectral
# IOP and depth and subtract the model value from observed Rrs using AM03.
#==================================================================================================

Saber_retrieve_rb_wise <-  function(use_true_IOPs = T,
                               a_non_water_path = "./data/rb_retrieve_demo_a.csv", 
                               bb_non_water_path = "./data/rb_retrieve_demo_bb.csv",
                              
                               chl=4.96, acdom440=0.9322314, anap440=0.07, 
                               bbp.550=0.00726002, 
                               
                               slope.parametric = FALSE, dg_composite = FALSE,
                               plot = FALSE, verbose = FALSE,
                               
                               z=2, 
                               
                               waves_sim = seq(400,800,10), obs_rrs
){
  
  ## OAC and IOP initialization
  
  # Read input parameters for forward modeling
  
  base.CDOM <- acdom440 # CDOM absorption at 440nm [m^-1]
  C_ph <- chl           # Phytoplankton concentration [mg/m^-3]
  base.NAP <- anap440   # NAP absorption at 440nm [m^-1]
  base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  base_CDM <- base.CDOM + base.NAP #CDM absorption at 440 nm [m^-1]
  lambda <- waves_sim       #Desired spectral range for simulation
  
  ##If in situ IOPs are provided
  if (!is.null(a_non_water_path) & !is.null(bb_non_water_path)) {
    
    a_non_water_df = read.csv(file = a_non_water_path, header = T)
    bb_non_water_df = read.csv(file = bb_non_water_path, header = T)
    
    if (is.data.frame(a_non_water_df) & is.data.frame(bb_non_water_df)) {
      
      a_non_water = a_non_water_df$at_w
      a_non_water_wave = a_non_water_df$wave
      
      bb_non_water = bb_non_water_df$bbp
      bb_non_water_wave = bb_non_water_df$wave
      
    } else {
      print("The IOPs must be provided as data frames with <<wave>> as wavelength & <<a/bb>> as value")
    }
    
  }
  
  #Shallow params
  zRb <- z             #bottom depth
  
  obs_rrs <- obs_rrs #Actual AOP data for which we try to retrieve Rb
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  ## Pure water absorption (1/m)
  # wavelength range [190;4000] [nm]
  abs.water <- read.table("./data/input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  if (verbose == TRUE) {
    print(paste0("Water absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm.")) 
  }
  
  ## Pure water backscattering (1/m)
  
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
  ## Absorption
  #--------------------------------------------------------------------------
  if (use_true_IOPs == FALSE) {
    print(paste0("Full spectral IOPs not provided, bio-optical models are used"))
    
    ## Plankthon absorption (1/m)
    # load plankton absorption data
    A0_A1_PhytoPlanc <- read.table("./data/input-spectra/A0_A1_PhytoPlanc.dat")
    # extract the values from the table
    lam_p <- A0_A1_PhytoPlanc$V1
    a0_p <- A0_A1_PhytoPlanc$V2
    a1_p <- A0_A1_PhytoPlanc$V3
    
    a0 <- rep(0,length(lambda))# [m^2/mg]
    a1 <- rep(0,length(lambda))# [m^2/mg]
    
    a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = lambda,method = "linear")$y
    a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = lambda,method = "linear")$y
    
    # Compute the value of plankton absorption as function of the concentration and wavelength
    aph_440 <- 0.06*(C_ph)^0.65# [mg/m^3]
    abs_ph <- rep(0,length(lambda))
    
    for (i in 1:length(lambda)){
      
      abs_ph[i] <- (a0[i] + a1[i]*log(aph_440))*aph_440
      
    }
    abs_ph[abs_ph < 0] <- 0
    
    if (verbose == TRUE) {
      print(paste0("Planktonic absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    #Calculate absorption of CDOM and NAP
    if (dg_composite == TRUE) {
      
      if (slope.parametric == TRUE) {
        if (type_Rrs_below_rb == "deep") {
          
          #parametric formula to retrieve spectral slope of CDOM + NAP
          S_CDM = 0.015 + (0.002/(0.6 + (obs_rrs[which.min(abs(lambda - 443))]/obs_rrs[which.min(abs(lambda - 555))])))
          
          print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
          
        } else {
          
          print("The water type is shallow, QAA will obtain wrong s_CDOM, thus const. 0.017 is used")
          S_CDM <- 0.017 #<< USER INPUT >>
          
        }
      } else {
        
        print("const. 0.017 is used as s_dg")
        S_CDM <- 0.017 #<< USER INPUT >>
      }
        #Calculate the absorption of dg
        Ga_CDM <- 1# [m^2/mg]
        Oa_CDM <- 0
        
        abs_CDM_440 <-  (Ga_CDM*base_CDM)+Oa_CDM# [1/m], CDOM abs. coeff. at 440 [nm]
        
        abs_CDM <- rep(0,length(lambda))
        
        for (i in 1:length(lambda)){
          
          abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
          
        }
    } else {
      ## CDOM absorption coefficient [1/m]
      
      Ga_CDOM <- 1# [m^2/mg]
      Oa_CDOM <- 0
      S_CDOM <- 0.014 #<< USER INPUT >>
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
      S_X <- 0.01160 #<< USER INPUT >>
      
      abs_X_440 <-  (Ga_X*base.NAP)+Oa_X# [1/m], SPM abs. coeff. at 440 [nm]
      
      abs_X <- rep(0,length(lambda))
      
      for (i in 1:length(lambda)){
        
        abs_X[i]  <- abs_X_440*exp(-S_X*(lambda[i] - 440))
        
      }
      if (verbose == TRUE) {
        print(paste0("NAP absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      abs_CDM = abs_CDOM + abs_X
    }
    
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
    
    if (plot == "TRUE") {
      ggsave(paste0("./gfx/forward_runs/SABER.absorption.model_chl=",signif(C_ph, digits = 3) ,"_acdom440=", signif(base.CDOM, digits = 3),"_anap440=", signif(base.NAP,digits = 3), ".png"), plot = g,
             scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
    }
    
    #--------------------------------------------------------------------------
    ## Scattering and back-scattering
    #--------------------------------------------------------------------------
    
    if (slope.parametric == TRUE) {
      #parametric formula to retrieve spectral slope of bbp
      refexponent = 2*(1-(1.2*exp(-0.9 * (obs_rrs[which.min(abs(lambda - 443))]/obs_rrs[which.min(abs(lambda - 555))]))))
      print(paste0("The spectral slope for bbp is calculated as:", refexponent))
      
    } else {
      refexponent <- 0.46 # << USER INPUT >>
      print(paste0("The spectral slope for bbp is constant as:", refexponent))
    }
    
    bb_x <- rep(0,length(lambda))# Backscattering coefficient for suspended particles [1/m]
    
    for (i in 1:length(lambda)){
      
      #print(lambda[i])
      bb_x[i] = base.bbp*((lambda[i]/550)^-(refexponent)) #Implement power model
      
    }
    if (verbose == TRUE) {
      print(paste0("NAP backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    #--------------------------------------------------------------------------
    ## Water system total IOPs
    #--------------------------------------------------------------------------
    
    ## Total Absorption Coefficient (1/m)
    
    abs <-  a_W + abs_ph + abs_CDM
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
        sensor_wl = c(394, 420, 470, 532, 620, 700)
        x = 555/sensor_wl
        print("HS-6 VSF is used")
      }
      if (length(bb_non_water_wave) == 9) {
        sensor_wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
        x = 555/sensor_wl
        print("BB-9 VSF is used")
      }
      
      #nz=length(IOP.fitted.down$Depth)
      nz=1 #As we only are interested for surface bb
      bbP555.down =rep(0,nz)
      nuP.down    =rep(0,nz)
      
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
  Ars1 <- 1.1576; Ars2 <- 1.0389 #Parametric coefficients for shallow water
  
  
  Rrs_Bottom <- (obs_rrs - Rrs_below_deep * (1 - Ars1 * exp(-zRb * (Kd + kuW)))) / (Ars2 * exp(-zRb * (Kd + kuB)))
  print(paste0("Bottom reflectance calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  print("!!!! The Bottom reflectance is not scaled with pi !!!!")
  return(data.frame("wavelength"=lambda, "rb"=Rrs_Bottom))
}

#=============================================================
## BATCH Rb retrieval for WISE-Man in situ
#=============================================================
Saber_generate_rb_batch <- function(staionlist, save_plot=T, save_data=T) {
  
  shallow_qc_stations = staionlist
  
  for (i in 1:length(shallow_qc_stations)) {
    stat_temp = shallow_qc_stations[i]
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;36m","##################### Rb retrieval start for Station: ", stat_temp, " #########################","\033[0m","\n"))
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    
    iop_aop_rb = suppressWarnings(read.surface.IOPs.wise(station.args = stat_temp,save_on_disc = T, 
                                                         verbose = F))
    
    bgc_params_rb = get_in_situ_params(station_name = stat_temp, use_bb_nup = T)
    
    IOP_files = list.files("./data/Rb_spectral/surface_iops/", full.names = T)
    
    idx_a_rb = grep(IOP_files, pattern = paste0("abs_surf_",stat_temp, ".csv$"))
    idx_bb_rb = grep(IOP_files, pattern = paste0("bb_surf_",stat_temp, ".csv$"))
    
    ###======================================================================
    # 4. Retrieve spectral Rb by providing IOPs and Depth
    ###======================================================================
    if (type_Rrs_below_rb == "shallow") {
      #Generate Rb using forward model parameters
      if (insitu.present == FALSE) {
        
        Rrs_bottom_gen = Saber_generate_rb(lambda = wavelength, fA = fA.set)
        print(paste0("Bottom composition is linear mixture of sand = ", fA.set[1], " , sediment = ",
                     fA.set[2], " , Chara contraria = ", fA.set[3],
                     " , Potamogeton perfoliatus = ", fA.set[4],
                     " , Potamogeton pectinatus = ", fA.set[5]))
        
        print("in situ Rrs is not supplied, user simulated Rrs is used for bottom reflectance retrieval")
        
        if (type_Rrs_water == "below_surface") {
          
          obs_rrs =  rrs.forward.am.param.conc.true_iop
        } else {
          obs_rrs = surface_rrs_translate(Rrs = rrs.forward.am.param.conc.true_iop)
        }
        
      } else {
        print("Historic data of the bottom composition is unknown for the station")
        print("in situ Rrs is supplied for bottom reflectance retrieval")
        obs_rrs = surface_rrs_translate(Rrs = iop_aop_rb$Rrs_0p)
      }
      
      cat(paste0("\033[0;36m","=======================PARAMETRIC RETIREVAL=================================================","\033[0m","\n"))
      #Retrieve spectral Rb from forward SABER generated Rrs with BGC parametrization
      Rrs_bottom_est_bgc = Saber_retrieve_rb_wise(use_true_IOPs = F, 
                                                  a_non_water_path = IOP_files[idx_a_rb],
                                                  bb_non_water_path = IOP_files[idx_bb_rb],
                                                  
                                                  chl = bgc_params_rb$chl_invivo, 
                                                  acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
                                                  anap440 = bgc_params_rb$abs_invivo$a_nap,
                                                  bbp.550 = bgc_params_rb$bbp555,
                                                  
                                                  slope.parametric = F, dg_composite = F,
                                                  
                                                  #z = 2,
                                                  z = iop_aop_rb$zB_COPS,
                                                  
                                                  #obs_rrs = rrs.forward.am.param.conc.true_iop,
                                                  obs_rrs =  obs_rrs
      )
      
      cat(paste0("\033[0;36m","=======================IOP based RETIREVAL=================================================","\033[0m","\n"))
      #Retrieve spectral Rb from forward SABER generated Rrs with True IOPs
      Rrs_bottom_est_iop = Saber_retrieve_rb_wise(use_true_IOPs = T, 
                                                  a_non_water_path = IOP_files[idx_a_rb],
                                                  bb_non_water_path = IOP_files[idx_bb_rb],
                                                  
                                                  chl = bgc_params_rb$chl_invivo, 
                                                  acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
                                                  anap440 = bgc_params_rb$abs_invivo$a_nap,
                                                  bbp.550 = bgc_params_rb$bbp555,
                                                  
                                                  slope.parametric = F, dg_composite = F,
                                                  
                                                  #z = 2,
                                                  z = iop_aop_rb$zB_COPS,
                                                  
                                                  #obs_rrs = rrs.forward.am.param.conc.true_iop,
                                                  obs_rrs = obs_rrs
      )
      
      cat(paste0("\033[0;36m","=======================SAVING OUTPUT=================================================","\033[0m","\n"))
      #Compare the actual vs retrieved Rb as a function of wavelength (#Change the plot to ggplot)
      if (insitu.present == TRUE) {
        bottom_ref_df = data.frame("wave"= as.numeric(wavelength),
                                   "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
                                   "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
                                   "rrs_shallow" = obs_rrs) 
        
        xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
        ymin= 0; ymax=1.25*max(bottom_ref_df$rb_est_iop[bottom_ref_df$wave < "650"], na.rm = T)
        ystp= signif(ymax/5, digits = 2)
        asp_rat <- (xmax-xmin)/(ymax-ymin)
        #browser()
        g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
                                   size=1.3,show.legend = TRUE) +
          geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
                    size=1.3,show.legend = TRUE)+
          geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
                    size=1.3,show.legend = TRUE)+
          
          scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
                                         expression(paste(italic("R")["B, model_BOPT"])),
                                         expression(paste(italic("R")["B, model_IOP"]))
          ), 
          values = c("purple","navyblue", "red")) +
          ggtitle(paste0(stat_temp)) +
          scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
                             limits = c(xmin, xmax), 
                             breaks = seq(xmin, xmax, xstp))  +
          scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
                             breaks = seq(ymin, ymax, ystp))+
          coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                      ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
          theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
                axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.ticks.length = unit(.25, "cm"),
                legend.position=c(0.6,1),
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
        
        
      } else {
        bottom_ref_df = data.frame("wave"= as.numeric(wavelength), 
                                   "rb_actual" = as.numeric(Rrs_bottom_gen$rb),
                                   "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
                                   "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
                                   "rrs_shallow" =obs_rrs) 
        
        xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
        ymin= 0; ymax=1.25*max(bottom_ref_df$rb_actual[bottom_ref_df$wave < "650"], na.rm = T)
        ystp= signif(ymax/5, digits = 2)
        asp_rat <- (xmax-xmin)/(ymax-ymin)
        
        g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
                                   size=1.3,show.legend = TRUE) +
          geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
                    size=1.3,show.legend = TRUE)+
          geom_point(data = bottom_ref_df,aes(y=rb_actual,x = wave,color="xx2"), 
                     size=3,show.legend = TRUE)+
          geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
                    size=1.3,show.legend = TRUE)+
          
          scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
                                         expression(paste(italic("R")["B, actual"])),
                                         expression(paste(italic("R")["B, model_BOPT"])),
                                         expression(paste(italic("R")["B, model_IOP"]))
          ), 
          values = c("purple","black","navyblue", "red")) +
          ggtitle(paste0(stat_temp)) +
          scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
                             limits = c(xmin, xmax), 
                             breaks = seq(xmin, xmax, xstp))  +
          scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
                             breaks = seq(ymin, ymax, ystp))+ 
          coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                      ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
          theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
                axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.ticks.length = unit(.25, "cm"),
                legend.position=c(0.6,1),
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
      }
      
      g 
      cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
      cat(paste0("\033[0;32m","##################### Rb retrieval end for Station: ", stat_temp, " #########################","\033[0m","\n"))
      cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
      
      #Save the bottom ref to disk
      if (save_data == TRUE) {
        write.csv(file = paste0("./Outputs/Bottom_ref/Rb_spectral_data_",stat_temp,".csv"), 
                  x = bottom_ref_df, quote = F, sep = ",", col.names = T, row.names = F)
      }
      
      
      if (save_plot == TRUE) {
        ggsave(filename = paste0("./Outputs/Bottom_ref/Rb_spectral_station_",stat_temp,".png"), plot = g,
               scale = 1.5, width = 4.5, height = 4.5, units = "in", dpi = 300)
      }
    } else {
      print("Rb retrieval is not applicable for deep water column")
    }
    
  }
}

