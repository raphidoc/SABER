#' rrs_bottom_lmm
#' Compute bottom remote sensing reflectance based on Linear Mixture Model
#'  of end member spectra.
#'
#'  @param rrs_end_member a tibble with columns {wavelength, class, rrs_b} end members spectrum
#'  @param rb_fraction a tibble with columns {class, fraction}, fraction must sum to 1
#'  @param wavelength requested interpolation wavelength
#'
#' @return a tibble with columns {wavelength, rrs_bottom}
#'  @export

rrs_bottom_lmm <- function(
    rrs_end_member,
    end_member_fraction,
    wavelength,
    verbose = F
    ) {

  # Reflection factors of bottom surface [1/sr]
  B1 <- 1/pi
  B2 <- 1/pi
  B3 <- 1/pi
  BOTTOM <- c(B1,B2,B3)

  # Bottom Albedo Calculation
  abott1 <-  rrs_end_member$class1
  abott2 <-  rrs_end_member$class2
  abott3 <-  rrs_end_member$class3

  abott <- rbind(abott1, abott2, abott3)

  areal_fraction <- end_member_fraction$fraction

  Bottom <-  matrix(nrow = length(areal_fraction), ncol = ncol(abott), 0)

  Bottom = areal_fraction * abott

  Rrs_Bottom = BOTTOM * Bottom # Bottom Rrs [1/sr]

  Bottom <- colSums(Bottom) # [unitless]


  rrs_bottom <- colSums(Rrs_Bottom) # [1/sr]

  if (!identical(rrs_end_member$wavelength, wavelength)) {
    rrs_bottom_interpolated <- approx(rrs_end_member$wavelength, rrs_bottom, wavelength)$y

    if (verbose) {
      plot(rrs_end_member$wavelength, rrs_bottom, xlab="wavelength",
           ylab="non-water absorption [m^-1]")
      lines(wavelength, rrs_bottom_interpolated, col="red", lwd=3)
    }

    rrs_bottom <- rrs_bottom_interpolated
  }

  return(
    tibble(
      wavelength,
      rrs_bottom
      )
    )
}
