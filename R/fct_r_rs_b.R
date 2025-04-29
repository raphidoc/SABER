#' Load bottom reflectance data
#' @param wavelength numeric vector
#' @param rrs_b matrix (columns = class names)
#' @export
load_r_rs_b <- function(wavelength, rrs_b) {
  stopifnot(is.numeric(wavelength), is.matrix(rrs_b))
  .Call("c_load_r_rs_b", as.numeric(wavelength), rrs_b)
}

#' Compute interpolated Rrs_b mixture
#' @param fractions named numeric vector (e.g., c(sand = 0.6, seagrass = 0.4))
#' @param wavelength numeric vector
#' @export
compute_r_rs_b_lmm <- function(fractions, wavelength) {
  .Call("c_compute_r_rs_b_lmm", fractions, as.numeric(wavelength))
}
