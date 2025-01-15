#' refraction_law
#' 
#' @param incidence_angle angle of incidence of electromagnetic radiation (rad)
#' @param n_1 refractive index of media of incidence (e.g. n_air = 1)
#' @param n_2 refractive index of media of transmittance (e.g. n_water = 1.34)
#' 
#' @author Soham Murkerjee
#' 
#' @references https://en.wikipedia.org/wiki/Snell%27s_law
#' 
#' @export

refraction_law <- function(incidence_angle, n1, n2){
  
  transmitance_angle = asin((n1 / n2) * sin(incidence_angle))
  
  return(transmitance_angle)
}