#' Fresnel equation for unpolarized light (Jerlov 1976)
#' 
#' @author Soham Mukherjee, Raphael Mabit
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

fresnel_reflectance <- function(incidence_angle, transmitance_angle){
  
  # Fresnel Law
  reflectance_factor = 0.5 * abs(
    (
      sin(incidence_angle-transmitance_angle)^2
      /
        sin(incidence_angle+transmitance_angle)^2
    )
    +
      (
        tan(incidence_angle-transmitance_angle)^2
        /
          tan(incidence_angle+transmitance_angle)^2
      )
  )
  
  return(reflectance_factor)
} 