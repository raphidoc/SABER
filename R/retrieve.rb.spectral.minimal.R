#==================================================================================================
# Saber_retrieve_rb_wise_minimal.R performs spectral retrieval of bottom reflectance (Rb) based on 
# Radiative transfer in optically complex case-II waters with miimal calculation for vectorization. 
# The Rb is retrieved from the Remote-Sensing Reflectance (Rrs) in optically complex waters from
# user-given spectral IOP and depth and subtract the model value scaled with attenuation coeffs.
# from observed Rrs using AM03.
#==================================================================================================

Saber_retrieve_rb_wise_minimal <-  function(
                                    obs_rrs,
                                    z,
                                    
                                    Rrs_deep,
                                    
                                    Kd, 
                                    KuW, 
                                    KuB,
                                    
                                     
                                    
                                    waves_sim = seq(400,800,10), 
                                    plot = FALSE, 
                                    verbose = FALSE
                                    
){
  
  ## OAC and IOP initialization
  
  # Read input parameters for forward modeling
  
  lambda <- waves_sim #Desired spectral range for simulation
  
  
  #Shallow params
  zRb <- z             #bottom depth
  
  obs_rrs <- obs_rrs #Actual AOP data for which we try to retrieve Rb
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #Final calculation for shallow Rrs
  Ars1 <- 1.1576; Ars2 <- 1.0389 #Parametric coefficients for shallow water
  
  
  Rrs_Bottom <- (obs_rrs - Rrs_below_deep * (1 - Ars1 * exp(-zRb * (Kd + kuW)))) / (Ars2 * exp(-zRb * (Kd + kuB)))
  
  if (verbose == T) {
    print(paste0("Bottom reflectance calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    print("!!!! The Bottom reflectance is not scaled with pi !!!!")
  }
  
  
  return(list("rb"=Rrs_Bottom))
}
