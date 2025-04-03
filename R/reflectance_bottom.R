#' rrs_bottom_lmm
#' Compute bottom remote sensing reflectance based on Linear Mixture Model
#'  of end member spectra.
#'
#'  @param rb a tibble with columns {wavelength, class, rrs_b} end members spectrum
#'  @param rb_fraction a tibble with columns {class, fraction}, fraction must sum to 1
#'  @param wavelength requested interpolation wavelength
#'
#' @return a tibble with columns {wavelength, rrs_bottom}
#'  @export

rrs_bottom_lmm <- function(
    rrs_end_member = read_csv(fs::path_package("SABER", "extdata", "rb_endmembers_algae-wise.csv")),
    end_member_fraction,
    wavelength
    ) {

  # Reflection factors of bottom surface [1/sr]
  B1 <- 1/pi
  B2 <- 1/pi
  B3 <- 1/pi
  BOTTOM <- c(B1,B2,B3)

  # Bottom Albedo Calculation
  abott1 <-  rb$class1
  abott2 <-  rb$class2
  abott3 <-  rb$class3

  abott <- rbind(abott1, abott2, abott3)

  Bottom <-  matrix(nrow = length(fA), ncol = ncol(abott), 0)
  Rrs_Bottom <- matrix(nrow = length(fA), ncol = ncol(abott), 0) # Bottom Rrs [1/sr]

  Bottom = fA * abott

  Rrs_Bottom = BOTTOM * Bottom

  Bottom <- colSums(Bottom) # [unitless]
  Rrs_Bottom <- colSums(Rrs_Bottom) # [1/sr]

  return(rrs_botom)
}
