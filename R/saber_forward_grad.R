
Saber_forward.grad <-  function(pars, bbp.550 = Fit.input$bbp.550
                           # plot = FALSE, verbose = FALSE,
                           # realdata = demo.rrs$rrs.demo, 
                           # realdata.exist = TRUE
                           ){
  
  ## OAC and IOP initialization
  
  # Read input params for forward modeling
    base.CDOM <- pars[2] # CDOM absorption at 440nm [m^-1]
    C_ph <- pars[1]           # Phytoplankton concentration [mg/m^-3]
    base.NAP <- pars[3]   # NAP absorption at 440nm [m^-1]
    base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
    lambda <- seq(400,800,10) 

   #Desired spectral range for simulation
  
  # #Shallow params
  # if (type_Rrs_below == "shallow") {
  #   base.CDOM <- pars[2] # CDOM absorption at 440nm [m^-1]
  #   C_ph <- pars[1]           # Phytoplankton concentration [mg/m^-3]
  #   base.NAP <- pars[3]   # NAP absorption at 440nm [m^-1]
  #   base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  #   lambda <- wavelength 
  #   fA <- fA.set     #Aerial fraction of bottom albedo
  #   #fA <- c(rb0,rb1,rb2,rb3,rb4,rb5)     #Aerial fraction of bottom albedo
  #   zB <- zB              #bottom depth
  #   
  # }
  
  #Rrs_obs.interp <- realdata #Actual AOP data for which we try to simulate the forward model
  #lambda <- wavelength
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Absorption
  #--------------------------------------------------------------------------
  ## Pure water absorption (1/m)
  
  # wavelength range [190;4000] [nm]
  abs.water <- read.table("./input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  
  ## Plankton absorption (1/m)
  
  # load plankton absorption data
  A0_A1_PhytoPlanc <- read.table("./input-spectra/A0_A1_PhytoPlanc.dat")
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
  #abs_ph[abs_ph < 0] <- 0
  
  
  ## CDOM absorption coefficient [1/m]
  
  Ga_CDOM <- 1# [m^2/mg]
  Oa_CDOM <- 0
  S_CDOM <- 0.014
  abs_CDOM_440 <-  (Ga_CDOM*base.CDOM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
  
  abs_CDOM <- rep(0,length(lambda))
  
  for (i in 1:length(lambda)){
    
    abs_CDOM[i]  <- abs_CDOM_440*exp(-S_CDOM*(lambda[i] - 440))
    
  }
  
  ## SPM absorption coefficient [1/m]
  
  Ga_X <- 1# [m^2/mg]
  Oa_X <- 0
  S_X <- 0.01160
  abs_X_440 <-  (Ga_X*base.NAP)+Oa_X# [1/m], SPM abs. coeff. at 440 [nm]
  
  abs_X <- rep(0,length(lambda))
  
  for (i in 1:length(lambda)){
    
    abs_X[i]  <- abs_X_440*exp(-S_X*(lambda[i] - 440))
    
  }
  
  
  #--------------------------------------------------------------------------
  ## Scattering and back-scattering
  #--------------------------------------------------------------------------
  
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
  
  
  # ## SPM backscattering coeff
  
  # bbxS_W <- 0.0086# [m^2/g] specfic backscattering according to WASI
  # bbxS <- ((bbxS_W)*(33.57*(10^(-6))))/g_size# [m^2/g] specfic backscattering as function of grain size
  # # b_ratio = 0.019;                                 # ratio between specific scattering and specific backscatteirng
  # # bxS = b_ratio*bbxS;                              # [m^2/g] specific scattering
  
  refexponent <- 0.46
  bb_x <- rep(0,length(lambda))# Backscattering coefficient for suspended particles [1/m]
  
  for (i in 1:length(lambda)){
    
    #print(lambda[i])
    bb_x[i] = base.bbp*((lambda[i]/550)^-(refexponent)) #Implement power model
    
  }
  #--------------------------------------------------------------------------
  ## Water system IOPs
  #--------------------------------------------------------------------------
  
  ## Total Absorption Coefficient (1/m)
  
  abs <-  a_W + abs_ph + abs_CDOM + abs_X
  #abs= a_W + a_tw;
  
  
  ## Total Backscattering Coefficient (1/m)
  
  bb <-  bb_W + bb_x
  #bb= bb_W + b_btw;
  
  ## Extinction Coefficient (1/m) and Single Back Scattering Albedo
  
  ext <-  abs+bb# [1/m] extinction coeff.
  omega_b <-  bb/ext# single back scattering albedo
  
  
  ################################################################################################
  ##                                      RT Model                                              ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance below the surface
  #--------------------------------------------------------------------------
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w; view_w <- geometry$view_w; rho_L <- geometry$rho_L
  
  
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
   } 
    
    if (type_Rrs_below == "shallow") {
      #Bottom Contribution

      # Reflection factors of bottom surface [1/sr]
      #B0 <- 1/pi;
      B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi; B4 <- 1/pi; B5 <- 1/pi;
      #BOTTOM <- c(B0,B1,B2,B3,B4,B5)
      BOTTOM <- c(B1,B2,B3,B4,B5)

      # Bottom Albedo (costant)
      # wavelength range [350;900] [nm]
      bott0<- read.table("./input-spectra/Bott0const.R")
      wavebottom <- bott0$V1
      Bott0 <-  bott0$V2
      abott0 <- rep(0,length(lambda))
      abott0 <-  Hmisc::approxExtrap(wavebottom, Bott0, xout = lambda, method = "linear")$y

      # Bottom Albedo Sand
      # wavelenght range [350;1000] [nm]
      bott1 <- read.table("./input-spectra/Bott1SAND.R")
      wavebottom <- bott1$V1
      Bott1 <-  bott1$V2
      abott1 <- rep(0,length(lambda))
      abott1 <-  Hmisc::approxExtrap(wavebottom, Bott1, xout = lambda, method = "linear")$y

      # Bottom Albedo of fine-grained sediment
      # wavelenght range [350;900] [nm]
      bott2 <- read.table("./input-spectra/Bott2silt.R")
      wavebottom <- bott2$V1
      Bott2 <-  bott2$V2
      abott2 <- rep(0,length(lambda))
      abott2 <-  Hmisc::approxExtrap(wavebottom, Bott2, xout = lambda, method = "linear")$y

      # Bottom Albedo of green makrophyte "Chara contraria"
      # wavelenght range [350;900] [nm]
      bott3 <- read.table("./input-spectra/Bott3chara.R")
      wavebottom <- bott3$V1
      Bott3 <-  bott3$V2
      abott3 <- rep(0,length(lambda))
      abott3 <-  Hmisc::approxExtrap(wavebottom, Bott3, xout = lambda, method = "linear")$y

      # Bottom Albedo of green makrophyte "Potamogeton perfoliatus"
      # wavelenght range [350;900] [nm]
      bott4 <- read.table("./input-spectra/Bott4perfol.R")
      wavebottom <- bott4$V1
      Bott4 <-  bott4$V2
      abott4 <- rep(0,length(lambda))
      abott4 <-  Hmisc::approxExtrap(wavebottom, Bott4, xout = lambda, method = "linear")$y

      # Bottom Albedo of green makrophyte "Potamogeton pectinatus"
      # wavelenght range [350;900] [nm]
      bott5 <- read.table("./input-spectra/Bott5pectin.R")
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



      #Final calculation for shallow Rrs
      Ars1 <- 1.1576; Ars2 <- 1.0389 #Paramteric coeffs for shallow water
      Rrs_below_shallow <-  Rrs_below_deep*(1-(Ars1*exp(-zB*(Kd+kuW)))) + Ars2*Rrs_Bottom*exp(-zB*(Kd+kuB))
      Rrs_below <- Rrs_below_shallow

    }
  Rrs <-  Rrs_below
  return(Rrs)
  }


