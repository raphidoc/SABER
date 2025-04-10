#' iop_from_oac
#' compute Inherent Optical Proprieties from Optically Active constituents.
#'
#' @param oac a tibble with column {chl [mg/m^-3], ag_440, bb_550}, optionally,
#'  provide manual ag and bb slope columns {ag_s_g, ag_s_d, bb_gamma}
#' @param rrs optionnal, a tibble with columns {wavelength, rrs} used for
#'  parametric retrieval of a_g and bb slope
#'
#' @importFrom assertthat has_name
#'
#' @export

iop_from_oac <- function(
    wavelength,
    oac,
    rrs = NULL,
    optically_shallow = NULL,
    verbose = F
    ) {

  # Phyto concentration to absorption  --------------------------------------

  # Plankton absorption (1/m)
  # load plankton absorption data
  A0_A1_PhytoPlanc <- read.table(fs::path_package("SABER", "extdata","a0_a1_phyto.csv"), sep = ",")

  # extract the values from the table
  lam_p <- A0_A1_PhytoPlanc$V1
  a0_p <- A0_A1_PhytoPlanc$V2
  a1_p <- A0_A1_PhytoPlanc$V3

  a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = wavelength, method = "linear")$y # [m^2/mg]
  a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = wavelength, method = "linear")$y# [m^2/mg]

  aph_440 <- 0.06 * (oac$chl)^0.65  # [mg/m^3] #Prieur & Satyendranath (1981)

  a_phy <- sapply(1:length(wavelength), function(i) (a0[i] + a1[i] * log(aph_440)) * aph_440)

  if (any(a_phy < 0)) {
    #rlang::warn("Some a_phy inferiro to 0")
    a_phy[a_phy < 0] <- 0
  }

  # CDOM absorption ---------------------------------------------------------

  ## CDOM+NAP absorption coefficient [1/m]

  Ga_CDOM <- 1
  Oa_CDOM <- 0

  abs_CDM_440 <- (Ga_CDOM*oac$ag_440)+Oa_CDOM # [1/m], CDOM abs. coeff. at 440 [nm]

  # 3 possibilities for CDOM slope: manual, parametric, default

  if (has_name(oac, c("ag_s_g", "ag_s_d"))) {
    cdom_slope = oac$ag_s_g + oac$ag_s_d

  } else if (isTRUE(optically_shallow) && tibble::is_tibble(rrs)) {
    cdom_slope = 0.015 + (
      0.002 / (0.6 + (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
                        rrs$rrs[which.min(abs(rrs$wavelength - 555))])))

  } else {
    cdom_slope <- 0.017
  }

  a_g <- sapply(1:length(wavelength), function(i) abs_CDM_440 * exp(-cdom_slope * (wavelength[i] - 440)))

  # Particulate backscattering ----------------------------------------------

  # 3 possibilities for bb slope: manual, parametric, default
  if (has_name(oac, "bb_gamma")) {
    bb_gamma <- oac$bb_gamma

  } else if (isTRUE(optically_shallow) && tibble::is_tibble(rrs)) {
    bb_gamma <- 2 * (1 - (1.2 * exp(
      -0.9 * (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
                rrs$rrs[which.min(abs(rrs$wavelength - 555))]))))

  } else {
    bb_gamma <- 0.46
  }

  bb_p <- sapply(1:length(wavelength), function(i) oac$bbp_550 * ((wavelength[i] / 550) ^ -bb_gamma))

  if (verbose) {
    if (verbose) {
      plot(wavelength, a_phy, xlab="wavelength",
           ylab="a_phy [m^-1]")

      plot(wavelength, a_g, xlab="wavelength",
           ylab="a_g [m^-1]")

      plot(wavelength, bb_p, xlab="wavelength",
           ylab="bb_p [m^-1]")
    }
  }

  return(
    tibble(
      wavelength,
      a_phy,
      a_g,
      bb_p,
      "a" = a_phy + a_g,
      "bb" = bb_p
    )
  )
}
