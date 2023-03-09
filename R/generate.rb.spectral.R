Saber_generate_rb <- function(lambda = seq(400,800,10), fA= fA.set){
  #Bottom Contribution
  
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
  return(data.frame("wavelength"=lambda, "rb"=Rrs_Bottom))
}

# Rrs_bottom_gen = Saber_generate_rb(lambda = wavelength, fA = fA.set)
# 
# Rrs_bottom_est = Saber_retrieve_rb(use_true_IOPs = F, chl = Fit.input$chl, 
#                                    acdom440 = Fit.input$acdom.440, 
#                                    anap440 = Fit.input$anap.440,
#                                    bbp.550 = Fit.input$bbp.550, z = zB, 
#                                    slope.parametric = F, obs_rrs = rrs.forward.am)
# plot(Rrs_bottom_gen$rb, Rrs_bottom_gen$rb)
