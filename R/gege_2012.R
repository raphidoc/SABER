#' Compute Rayleigh scattering at the air-sea interface (sky glint)
#' 
#' @author Soham Mukerjee, Raphael Mabit
#' 
#' @param alpha angstrom exponent
#' @param p atmospheric pressure [mb]
#' @param RH relative humidity
#' @param Hoz scale height for ozone ? [cm] or [cm-atm]
#' @param WV scale height water vapor ? [cm] or [cm-atm]
#' 
#' @references
#' Gordon, H.R., Brown, J.W. and Evans, R.H. (1988)
#'  ‘Exact Rayleigh scattering calculations for use with the Nimbus-7 
#'  Coastal Zone Color Scanner’, Applied Optics, 27(5), p. 862. 
#'  Available at: https://doi.org/10.1364/AO.27.000862.
#'  
#'  Gege, P. (2012) 
#'  ‘Analytic model for the direct and diffuse components of downwelling 
#'  spectral irradiance in water’, Applied Optics, 51(9), p. 1407. 
#'  Available at: https://doi.org/10.1364/AO.51.001407.
#'  
#' @export
gege_2012 <- function(
  wavelength,
  theta_s, #sun angle of sensor above water surface
  # Atmospheric conditions
  # Irradiance intensities [1/sr]
  g_dd=0.02,
  g_dsr=1/pi,
  g_dsa=1/pi,
  # Intensities of light sources 
  f_dd= 1,
  f_ds= 1,
  # Air Mass
  AM = 1,
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
    x = SABER::F0$wavelength,
    y = SABER::F0$F0,
    xout = wavelength,
    method = "linear")$y
  
  # Oxygen absorption [1/cm]
  a_O2_res <- Hmisc::approxExtrap(
    x = SABER::a_O2$wavelength,
    y = SABER::a_O2$a_O2,
    xout = wavelength,
    method = "linear")$y
  
  # Ozone absorption [1/cm]
  a_O3_res <- Hmisc::approxExtrap(
    x = SABER::a_O3$wavelength,
    y = SABER::a_O3$a_O3,
    xout = wavelength,
    method = "linear")$y
  
  # Water vapour absorption [1/cm]
  a_WV_res <- Hmisc::approxExtrap(
    x = SABER::a_WV$wavelength,
    y = SABER::a_WV$a_WV,
    xout = wavelength,
    method = "linear")$y
  
  theta_s <- theta_s*(pi/180)
  
  theta_s_t <- refraction_law(theta_s, n1 = 1, n2 = 1.34) 
  fresnel_rho_sun = fresnel_reflectance(theta_s, theta_s_t)
  
  # Atmospheric path length -------------------------------------------------
  M <- 1/(cos(theta_s)+(0.50572*((90 + 6.079975-theta_s)^(-1.253))))
  # M corrected for non standard atmopsheric pressure
  M1 <-  (M*P)/1013.25
  # path length for ozone
  Moz <-  1.0035/ (((cos(theta_s)^2)+0.007)^0.5)
  
  # Air mass type
  # TODO: type_case_water is defined in main.R as an integer outside the
  #   function. How was this working with this scope ?
  # type_case_water <- 1
  # if (type_case_water == "1") { 
  #   AM <- 1
  # } else{
  #   AM <- 1
  # }
  
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
  Fa <-  1-(0.5*exp((B1+(B2*cos(theta_s)))*cos(theta_s)))
  
  # Downwelling irradiances -------------------------------------------------
  # Direct 
  Edd <-  F0_res*Tr*Taa*Tas*Toz*To*Twv*cos(theta_s)
  # Diffuse Rayleigh
  Edsr <- (1/2)*F0_res*(1-(Tr^(0.95)))*Taa*Tas*Toz*To*Twv*cos(theta_s)
  # Diffuse aerosol
  Edsa <-  F0_res*(Tr^(1/2))*Taa*(1-Tas)*Toz*To*Twv*cos(theta_s)*Fa
  # Diffuse downwelling irradiance above the surface
  Eds <-  Edsr+Edsa
  
  # Downwelling irradiance just above the water surface
  Ed_surface <-  (f_dd*Edd) + (f_ds*Eds)
  
# WASI sky radiances ------------------------------------------------------
# see: Gege, P. (2015) ‘The Water Colour Simulator WASI’.

  # Sky Radiance [mW/sr m^2 nm]
  Ls <-  (g_dd*Edd) + (g_dsr*Edsr) + (g_dsa*Edsa)
  
  # Sun specular reflected radiance
  L_sun_glint <- fresnel_rho_sun*g_dd*Edd
  
  # Sky specular reflected radiance
  L_sky_glint <- fresnel_rho_sun*(g_dsa*Edsa+g_dsr*Edsr)
  
  result <- data.frame(
    "L_sky" = Ls,
    "Ed_surface" = Ed_surface,
    "L_sun_glint" = L_sun_glint,
    "L_sky_glint" = L_sky_glint
  )
    
  return(result)
}