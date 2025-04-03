#=====================================================================================================
# Saber_forward_fast.R simulates the remote sensing reflectance given the wavelengths,
# the water components along with bathymetry and bottom reflectance (for shallow water).

# This code has all the fastest modes of bio-optical parametrih_wations for IOPs. For All B-OPT models,
# refer to Saber.forward.final().

# The SA algorithm is excluded for inelastic scattering equivalent Rrs. Refer to Saber.forward.final()
# for the inelastic scattering support.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=====================================================================================================
#SABER FORWARD Model FAST version


#' Lee 1998 forward model
#' Compute Rrs given initial parameters
#' Question For Soham:
#' Where does the parametric formula to retrieve spectral slope of CDOM + NAP comes from, QAA ?
#' Where does the a_w values comes from ?
#' Why is pure water backscattering dependent on water type ?
#'
#' @author Soham Mukherjee
#'
#' @param iop a tibble with column {wavelength, a, bb}
#' @param oac a tibble with column {chl [mg/m^-3], ag_440, bb_550}, optionally,
#'  to provide manual ag and bb slope columns {ag_s_g, ag_s_d, bb_gamma}
#' @param rb optional, a tibble with columns {wavelength, class`X`} end member
#'  reflectance spectra.
#' @param rrs optionnal, a tibble with columns {wavelength, rrs} used for
#' parametric retrieval of a_g and bb slope
#'
#' @importFrom Hmisc approxExtrap
#' @importFrom assertthat has_name
#' @references Lee, z. et al. (1998) ‘Hyperspectral remote sensing for shallow waters. I. A semianalytical model’, Applied Optics, 37(27), pp. 6329–6338. Available at: https://doi.org/10.1364/AO.37.006329.
#'

lee98_forward <-  function(
    iop = NULL,
    oac = NULL,
    water_type = 2,
    optically_shallow = T,
    h_w = 2, # bottom depth
    rb_fraction, # areal fraction of bottom types
    rb = read_csv(fs::path_package("SABER", "extdata", "rb_endmembers_algae-wise.csv")),
    view = 40,
    sun = 20,
    wavelength = seq(400,800,10),
    rrs = NULL,
    verbose = T
){

# Pure water IOPs ---------------------------------------------------------
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

# Define IOPs  ------------------------------------------------------------

  # If spectral IOP are provided use them as is
  # Else if only OAC are provided use bio-optical to derive spectral IOP
  if (tibble::is_tibble(iop)) {

    if (!identical(iop$wavelength, wavelength)) {
      rlang::abort(
        "iop walength different than requested simulation wavelength,
        use `parse_iop` to prepare iop data frame")

    } else {
      a_non_water <- iop$a
      bb_non_water <- iop$bb

      a_t <-  a_w + a_non_water
      bb_t <-  bb_w + bb_non_water

    }

  } else if (tibble::is_tibble(oac)) {

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
      rlang::warn("a_phy inferiro to 0")
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

    } else if (!optically_shallow && tibble::is_tibble(rrs)) {

      cdom_slope = 0.015 + (
        0.002 / (0.6 + (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
                          rrs$rrs[which.min(abs(rrs$wavelength - 555))])))

    } else {
      cdom_slope <- 0.017
    }

    #Vectorih_wation
    a_g <- sapply(1:length(wavelength), function(i) abs_CDM_440 * exp(-cdom_slope * (wavelength[i] - 440)))

# Particulate backscattering ----------------------------------------------

    # 3 possibilities for bb slope: manual, parametric, default
    if (has_name(oac, "bb_gamma")) {
      bb_gamma <- oac$bb_gamma
    } else if (!optically_shallow && tibble::is_tibble(rrs)) {
      bb_gamma <- 2 * (1 - (1.2 * exp(
        -0.9 * (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
                  rrs$rrs[which.min(abs(rrs$wavelength - 555))]))))
    } else {
      bb_gamma <- 0.46
    }

    bb_p <- sapply(1:length(wavelength), function(i) oac$bbp_550 * ((wavelength[i] / 550) ^ -bb_gamma))

# Total IOPS --------------------------------------------------------------

    ## Total Absorption Coefficient (1/m)
    a_t <-  a_w + a_phy + a_g

    ## Total Backscattering Coefficient (1/m)
    bb_t <-  bb_w + bb_p

  } else {
    rlang::abort("Please provide IOP or OAC")
  }

  ## Extinction Coefficient (1/m) and Single Back Scattering Albedo

  ext <-  a_t + bb_t # [1/m] extinction coeff.
  omega_b <-  bb_t / ext # single back scattering albedo

# Radiative Transfer Model ------------------------------------------------

  ## Remote sensing reflectance 0m the surface
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w
  view_w <- geometry$view_w
  rho_L <- geometry$rho_L

  ## Remote Sensing Reflectance 0m the water surface
  p1  <- 0.084
  p2  <- 0.17
  k1w <- 1.03
  k2w <- 2.04
  k1b <- 1.04
  k2b <- 5.04
  q1 <- 1 #(for viewing angle=0)
  k <- a_t + bb_t
  u <- bb_t / k

  rrs_0m_deep <- q1 * (p1 + p2 * u) * u

  if (optically_shallow) {

    #Shallow parameters
    fA <- rb_fraction     #Aerial fraction of bottom albedo

    # Reflection factors of bottom surface [1/sr]
    B1 <- 1/pi
    B2 <- 1/pi
    B3 <- 1/pi
    BOTTOM <- c(B1, B2, B3)#,B4,B5)

    # Bottom Albedo Calculation
    abott1 <-  rb$class1
    abott2 <-  rb$class2
    abott3 <-  rb$class3

    abott <- rbind(abott1, abott2, abott3)#, abott4, abott5)

    Bottom <-  matrix(nrow = length(fA), ncol = ncol(abott), 0)
    Rrs_Bottom <- matrix(nrow = length(fA), ncol = ncol(abott), 0)

    Bottom = fA * abott #BitWISE operation

    Rrs_Bottom = BOTTOM * Bottom #BitWISE operation # Bottom Rrs [1/sr]

    Bottom <- colSums(Bottom) #[unitless]
    Rrs_Bottom <- colSums(Rrs_Bottom)# [1/sr]

    # Attenuation Coefficients
    mu_s <- cos(sun_w)
    mu_v <- cos(view_w)
    du_w <- k1w * sqrt(1 + k2w * u)
    du_b <- k1b * sqrt(1 + k2b * u)

    rrs_0m_shallow <- rrs_0m_deep *
      (1 - exp(-(1 / mu_s + du_w / mu_v) * k * h_w)) +
      Rrs_Bottom * exp(-(1 / mu_s + du_b / mu_v) * k * h_w)

    rrs_0m <- rrs_0m_shallow

  } else {

    rrs_0m <- rrs_0m_deep

  }

  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------

  if (verbose) {
    plot(wavelength, rrs_0m, xlab="wavelength",
         ylab="Rrs 0m [m^-1]")
  }

  return(tibble(
    wavelength,
    rrs_0m
  ))

}
