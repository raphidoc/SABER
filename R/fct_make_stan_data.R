#' Make base data needed by STAN model
#'
#' @param ... additional fixed data expected by the model
#'
#' @export
make_stan_data_base <- function(
    wavelength,
    water_type, theta_view, shallow,
    bottom_class_names,
    ...
) {

  dots <- list(...)

  # ---- required extra parameters ----
  required <- c("a_gnap_s", "bb_p_gamma")

  missing <- setdiff(required, names(dots))
  if (length(missing) > 0) {
    stop(
      "Missing required parameters in `...`: ",
      paste(missing, collapse = ", ")
    )
  }

  a_gnap_s    <- as.numeric(dots$a_gnap_s)
  bb_p_gamma  <- as.numeric(dots$bb_p_gamma)

  luts <- make_luts(
    wavelength = wavelength
  )

  # --- enforce Stan model expectation: exactly 3 bottom classes ---
  if (length(bottom_class_names) != 3) {
    stop("This Stan model expects exactly 3 bottom_class_names (for simplex[3] and bottom_class_ids[3]).")
  }

  bottom_class_ids <- match(bottom_class_names, luts$classes)
  if (any(is.na(bottom_class_ids))) {
    missing <- bottom_class_names[is.na(bottom_class_ids)]
    stop("Bottom class not found in library: ", paste(missing, collapse = ", "))
  }
  bottom_class_ids <- as.integer(bottom_class_ids)

  # --- build basis Phi on the wavelength grid ---
  wl <- as.numeric(wavelength)
  n_wl <- length(wl)

  list(
    n_wl = n_wl,
    wavelength = wl,

    a_w  = luts$a_w,
    a0   = luts$a0,
    a1   = luts$a1,
    bb_w = luts$bb_w,

    n_class = length(luts$classes),
    r_b_mu_lib = luts$r_b_mu_lib,
    r_b_sd_lib = luts$r_b_sd_lib,
    bottom_class_ids = bottom_class_ids,

    water_type = as.integer(water_type),
    theta_view = as.numeric(theta_view),
    shallow = as.integer(shallow),

    a_gnap_s = a_gnap_s,
    bb_p_gamma = bb_p_gamma
  )
}

#' Make per observation data needed by STAN model
#'
#' @export
make_stan_data_obs <- function(
    df,
    stan_data_base,
    use_measured_sigma = TRUE,
    sigma_fallback = 0.00005,
    sigma_floor = 1e-8
    ) {

  wl <- as.numeric(df$wavelength)
  rrs_obs <- as.numeric(df$rrs_0m)
  theta_sun <- abs(unique(as.numeric(df$theta_sun)))
  h_w <- abs(unique(as.numeric(df$h_w)))

  stopifnot(length(wl) == stan_data_base$n_wl)
  stopifnot(all(abs(wl - stan_data_base$wavelength) < 1e-8))

  sigma_vec <- if (use_measured_sigma) as.numeric(df$rrs_unc) else rep(sigma_fallback, stan_data_base$n_wl)

  # cleanup
  sigma_vec[!is.finite(sigma_vec)] <- NA_real_
  if (anyNA(sigma_vec)) {
    s0 <- median(sigma_vec, na.rm = TRUE)
    sigma_vec[is.na(sigma_vec)] <- if (is.finite(s0)) s0 else sigma_fallback
  }
  sigma_vec <- pmax(sigma_vec, sigma_floor)

  list(
    rrs_obs = rrs_obs,
    sigma_rrs = sigma_vec,
    theta_sun = theta_sun,
    h_w = h_w,
    wl = wl
  )
}

make_stan_data <- function(stan_data_base, stan_data_obs) {
  append(stan_data_base, stan_data_obs)
}
