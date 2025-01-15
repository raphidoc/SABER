#' Transmit Rrs(0-) to Rrs(0+)
#' 
#' @references 
#' Gege, P. (2012)
#'  ‘Analytic model for the direct and diffuse components of downwelling
#'   spectral irradiance in water’, Applied Optics, 51(9), p. 1407.
#'    Available at: https://doi.org/10.1364/AO.51.001407.
#' 
#' @author Soham Mukherjee
#' 
#' @export

rrs_subsurface_to_surface <- function(
    rrs_subsurface,
    wavelength,
    # Sun-Sensor Geometry above surface
    view_above = 15, #viewing angle of sensor above water surface
    sun_above = 30, #sun angle of sensor above water surfac
    # Water surface roughness
    q_surf = pi #[sr]
) {
  
  # angles from deg to rad
  geometry_above = snell_law(view = view_above, sun = sun_above)
  
  view_rad=geometry_above$view_w # rad
  sun_rad=geometry_above$sun_w   # rad
  rho_L_above = geometry_above$rho_L #fresnel reflectance

  sigma <- 0.03
  # Refraction index of water
  nW <- 1.34
  # Transmission factor water > air
  rho_U <- 0.54
  # Q coefficient, surface roughness
  Q <- q_surf
  
  # Sub-surface Rrs (0-) is translated to surface (0+)
  rrs_surf = (((1-sigma)*(1-rho_L_above)/(nW^2))* (rrs_subsurface/(1-rho_U*Q*rrs_subsurface)))
  

}
