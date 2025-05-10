#' debug function

build_cache <- function(wavelength) {
  .Call("c_build_cache", as.numeric(wavelength))
}

get_a_w <- function() {
  .Call("c_get_a_w")
}
