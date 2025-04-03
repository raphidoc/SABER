#' pure_water_iop
#' compute pure water iop at requested wavelength
#'
#' @param wavelength vector of requested wavelength
#' @param water_type vector of requested wavelength
#'
#' @export

pure_water_iop <- function(
    wavelength,
    water_type = 2
) {
  # absorption (1/m)
  a_water <- read.table(fs::path_package("SABER", "extdata", "aw.csv"), header = F, sep = ",")
  a_w <-  approx(a_water$V1, a_water$V2, wavelength)$y # abs. of pure water [1/m]

  ## backscattering (1/m)
  if (water_type == 1) {
    b1 <-  0.00144#  [1/m]

  } else if (water_type == 2) {
    b1 <-  0.00111#  [1/m]
  } else {
    rlang::abort("Water type are limited to case 1 and 2")
  }
  lambda1 <-  500# [nm]
  bb_w <- b1*(wavelength/lambda1)^(-4.32)# [1/m]

  return(
    tibble(
      wavelength,
      a_w,
      bb_w
    )
  )
}
