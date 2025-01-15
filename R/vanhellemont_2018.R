#' calculate the sky reflectance (sky glint, rho_sky) at the air-water interface
#' @author Soham Mukherjee
#' 
#' @references
#'  Gege, P. (2012)
#'  ‘Analytic model for the direct and diffuse components of downwelling
#'   spectral irradiance in water’, Applied Optics, 51(9), p. 1407.
#'    Available at: https://doi.org/10.1364/AO.51.001407.
#'    
#'  Vanhellemont, Q. and Ruddick, K. (2018)
#'   ‘Atmospheric correction of metre-scale optical satellite data for
#'    inland and coastal water applications’, Remote Sensing of Environment,
#'     216, pp. 586–597.
#'      Available at: https://doi.org/10.1016/j.rse.2018.07.015.
#' 
#' @export

calc_rho_sky <- function(
    rrs_subsurface,
    wavelength,
    # Sun-Sensor Geometry above surface
    view_above = 15, #viewing angle of sensor above water surface
    sun_above = 30, #sun angle of sensor above water surfac
    # Water surface roughness
    q_surf = pi, #[sr]
    # Atmospheric conditions
    # Irradiance intensities [1/sr]
    g_dd=0.05,
    g_dsr=0,
    g_dsa=0,
    # Intensities of light sources 
    f_dd= 1,
    f_ds= 1,
    # Angstrom exponent
    alpha = 1.317,
    # Atmospheric pressure 
    P = 1013.25, # [mbar]
    # Relative Humidity
    RH = 0.60,
    # Scale height for ozone
    Hoz = 0.300, # [cm]
    # Scale height of the precipitate water in the atmosphere
    WV= 2.500 # [cm])
) {
  
  # Extraterrestrial solar irradiance [mW/m^2 nm]
  F0_res <- Hmisc::approxExtrap(
    x =F0$wavelength,
    y =F0$F0,
    xout = wavelength,
    method = "linear")$y
  
  # Oxygen absorption [1/cm]
  a_O2_res <- Hmisc::approxExtrap(
    x =a_O2$wavelength,
    y =a_O2$a_O2,
    xout = wavelength,
    method = "linear")$y
  
  # Ozone absorption [1/cm]
  a_O3_res <- Hmisc::approxExtrap(
    x = a_O3$wavelength,
    y = a_O3$a_O3,
    xout = wavelength,
    method = "linear")$y
  
  # Water vapour absorption [1/cm]
  a_WV_res <- Hmisc::approxExtrap(
    x =a_WV$wavelength,
    y =a_WV$a_WV,
    xout = wavelength,
    method = "linear")$y
  
  # Convert deg -> rad
  view_above <- view_above*(pi/180)
  sun_above <- sun_above*(pi/180)
  
  # Check consistency with codebase function
  # test <- snell_law(view_above, sun_above)
  
  view_under <- refraction_law(view_above, n_1 = 1, n_2 = 1.34) 
  sun_under <- refraction_law(sun_above, n_1 = 1, n_2 = 1.34) 
  rho_L_above = fresnel_reflection(view_above, view_under)
  
  # Downwelling Irradiance [mW/m^2 nm]
  
  
  # Atmospheric path length -------------------------------------------------
  M <- 1/(cos(sun_under)+(0.50572*((90 + 6.079975-sun_under)^(-1.253))))
  # M corrected for non standard atmopsheric pressure
  M1 <-  (M*P)/1013.25
  # path length for ozone
  Moz <-  1.0035/ (((cos(sun_under)^2)+0.007)^0.5)
  
  # Air mass type
  # TODO: type_case_water is defined in main.R as an integer outside the
  #   function. How was this working with this scope ?
  type_case_water <- 1
  if (type_case_water == "1") { 
    AM <- 1
  } else{
    AM <- 1
  }
  
  # single scattering albedo
  omega_a <- ((-0.0032*AM) + 0.972)*exp(RH*3.06*(10^-4))
  
  # aerosol scale height [km]
  Ha <- 1 
  # horizontal visibility [km]
  V <- 15
  # 
  beta <- 3.91*(Ha/V)
  tau_a <-  beta*((wavelength/550)^(-alpha))
  
  
  # Atmospheric transmittance spectra ---------------------------------------
  # Rayleigh scattering
  Tr <-   exp(-M1/((115.6406*(wavelength^(4)))-(1.335*(wavelength^(2)))))
  # Aerosol absorption
  Taa <-  exp(-(1-omega_a)*tau_a*M)
  # Aerosol scattering
  Tas <-  exp(-omega_a*tau_a*M)
  # Ozone absorption
  Toz <-  exp(-a_O3_res*Hoz*Moz)
  # Oxygen absorption
  To <-   exp((-1.41*a_O2_res*M1)/((1+(118.3*a_O2_res*M1))^0.45))
  # Water vapor absorption
  Twv <-  exp((-0.2385*a_WV_res*WV*M)/((1+(20.07*a_WV_res*WV*M))^0.45))
  
  
  # Aerosol forward scattering probability ----------------------------------
  # asymmetry factor of the aerosol scattering phase function
  if (alpha < 0) {
    assymetry_factor <- 0.82
  } else if (alpha > 1.2) {
    assymetry_factor <- 0.65
  } else {
    assymetry_factor <- -0.1417*alpha+0.82
  }
  # Empirical coefficient for the Aerosol forward scattering probability
  B3 <- log(1 - assymetry_factor)
  B2 <-  B3*(0.0783 +(B3*(-0.3824-0.5874*B3)))
  B1 <-  B3*(1.459 +(B3*(0.1595+0.4129*B3)))
  # Aerosol forward scattering probability
  Fa <-  1-(0.5*exp((B1+(B2*cos(sun_above)))*cos(sun_above)))
  
  
  # Downwelling irradiances -------------------------------------------------
  Edd <-  F0_res*Tr*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsr <- (1/2)*F0_res*(1-(Tr^(0.95)))*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsa <-  F0_res*(Tr^(1/2))*Taa*(1-Tas)*Toz*To*Twv*cos(sun_rad)*Fa
  
  # diffuse downwelling irradiance (sum of Rayleigh and aerosol)
  Eds <-  Edsr+Edsa
  
  # downwelling irradiance just above the water surface
  Ed <-  (f_dd*Edd) + (f_ds*Eds)
  
  # Sky Radiance [mW/sr m^2 nm]
  Ls <-  (g_dd*Edd) + (g_dsr*Edsr) + (g_dsa*Edsa)
  
  # Remote sensing reflectance above the surface [1/sr]
  # Sky Glint equivalent Rrs calculated
  # Rrs_above is in fact Rrs_sky_glint ?
  # What about Rrs_sun_glint ?
  Rrs_glint <-  rho_L_above*(Ls/Ed)
}

