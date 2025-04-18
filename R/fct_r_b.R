#' Compute bottom reflectance based on a Linear Mixture Model
#' of end member spectra.
#'
#' @param r_b_spectrum A tibble with columns {wavelength, class, r_b} for end member spectra
#' @param r_b_fraction A vector named by class of r_b_spectrum and value corresponding to fraction (sum to 1).
#' @param wavelength Target wavelengths for interpolation
#' @param verbose If TRUE, plots original and interpolated spectra
#'
#' @return A tibble with columns {wavelength, r_b}
#' @export
r_b_lmm <- function(
    r_b_spectrum,
    r_b_fraction,
    wavelength,
    verbose = FALSE
) {
  stopifnot(all(c("wavelength", "class", "r_b") %in% names(r_b_spectrum)))
  stopifnot(is.numeric(r_b_fraction), !is.null(names(r_b_fraction)))

  total <- sum(r_b_fraction)
  if (abs(total - 1) > 1e-6) {
    warning("Fractions do not sum to 1. They will be normalized.")
    r_b_fraction <- r_b_fraction / total
  }

  # Merge spectra with class fraction
  fraction_df <- tibble::tibble(
    class = names(r_b_fraction),
    fraction = as.numeric(r_b_fraction)
  )

  r_b_combined <- dplyr::inner_join(r_b_spectrum, fraction_df, by = "class")

  if (nrow(r_b_combined) == 0) {
    stop("No overlapping classes between spectra and fractions.")
  }

  # Compute weighted sum
  r_b <- r_b_combined %>%
    dplyr::mutate(weighted_r_b = (1 / pi) * r_b * fraction) %>%
    dplyr::group_by(wavelength) %>%
    dplyr::summarise(r_b = sum(weighted_r_b))

  # Interpolate if wavelength mismatch
  if (!all(wavelength %in% r_b$wavelength)) {
    r_b <- tibble(
      wavelength = wavelength,
      r_b = approx(
        x = r_b$wavelength,
        y = r_b$r_b,
        xout = wavelength)$y
    )
  }

  if (verbose) {
    ply <- plot_ly() %>%
      add_lines(
        data = r_b_combined,
        x = ~wavelength,
        y = ~r_b,
        color = ~class
      ) %>%
      add_lines(
        data = r_b,
        x = ~wavelength,
        y = ~r_b,
        name = "R_b_lmm"
      )
    print(ply)
  }

  return(r_b)
}
