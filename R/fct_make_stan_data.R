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
  required <- c("a_nap_star", "bb_p_star", "a_g_s", "a_nap_s", "bb_p_gamma")

  missing <- setdiff(required, names(dots))
  if (length(missing) > 0) {
    stop(
      "Missing required parameters in `...`: ",
      paste(missing, collapse = ", ")
    )
  }

  a_nap_star <- as.numeric(dots$a_nap_star)
  bb_p_star <- as.numeric(dots$bb_p_star)

  a_g_s <- as.numeric(dots$a_g_s)
  a_nap_s <- as.numeric(dots$a_nap_s)
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

    a_nap_star = a_nap_star,
    bb_p_star = bb_p_star,

    a_g_s = a_g_s,
    a_nap_s = a_nap_s,
    bb_p_gamma = bb_p_gamma
  )
}

# #' Make base data needed by STAN model (low-rank PCA mixture bottom prior)
# #'
# #' @export
# make_stan_data_base_gmm <- function(
#     wavelength,
#     water_type, theta_view, shallow,
#     bottom_class_names,
#     rb_obs_path = "/home/raphael/R/SABER/inst/extdata/r_b_gamache_obs.csv",
#     q = 5,
#     use_empirical_pi = F,
#     ...
# ) {
#   dots <- list(...)
#
#   luts <- make_luts_gmm(
#     wavelength = wavelength,
#     rb_obs_path = rb_obs_path,
#     q = q
#   )
#
#   # needed <- c("classes","K","q","rb_U","rb_mu_k","rb_sigma2_k","rb_pi",
#   #             "a_w","a0","a1","bb_w")
#   # miss <- setdiff(needed, names(luts))
#   # if (length(miss) > 0) stop("make_luts_gmm() missing fields: ", paste(miss, collapse = ", "))
#
#   if (length(bottom_class_names) < 1) stop("bottom_class_names must have length >= 1")
#   bad <- setdiff(bottom_class_names, luts$classes)
#   if (length(bad) > 0) stop("Bottom class not found in library: ", paste(bad, collapse = ", "))
#
#   # subset in requested order
#   rb_mu_k <- lapply(bottom_class_names, function(cl) luts$rb_mu_k[[cl]])         # each length n_wl
#   rb_tau_k <- lapply(bottom_class_names, function(cl) luts$rb_tau_k[[cl]])       # each length q
#   rb_L_k <- lapply(bottom_class_names, function(cl) luts$rb_L_k[[cl]])
#   sigma2_full <- setNames(luts$rb_sigma2_k, luts$classes)
#   rb_sigma2_k <- unname(sigma2_full[bottom_class_names])
#
#   # weights (empirical proportions from luts, or uniform)
#   K <- length(bottom_class_names)
#   if (use_empirical_pi) {
#     pi_full <- setNames(luts$rb_pi, luts$classes)
#     rb_pi <- unname(pi_full[bottom_class_names])
#     rb_pi <- rb_pi / sum(rb_pi)
#   } else {
#     rb_pi <- rep(1 / K, K)
#   }
#
#   wl <- as.numeric(wavelength)
#   n_wl <- length(wl)
#
#   # hyperparameters: fail fast if missing
#   hyp_needed <- c("a_nap_star_mu","a_nap_star_sd","bb_p_star_mu","bb_p_star_sd",
#                   "a_g_s_mu","a_g_s_sd","a_nap_s_mu","a_nap_s_sd",
#                   "bb_p_gamma_mu","bb_p_gamma_sd","h_w_mu","h_w_sd")
#   hyp_miss <- setdiff(hyp_needed, names(dots))
#   if (length(hyp_miss) > 0) stop("Missing required hyperparameters in `...`: ", paste(hyp_miss, collapse = ", "))
#
#   list(
#     n_wl = n_wl,
#     wavelength = wl,
#
#     a_w  = luts$a_w,
#     a0   = luts$a0,
#     a1   = luts$a1,
#     bb_w = luts$bb_w,
#
#     # low-rank mixture prior inputs
#     K = as.integer(K),
#     q = as.integer(luts$q),
#     rb_U = luts$rb_U,                                # matrix [n_wl, q]
#     rb_mu_k = rb_mu_k,                               # list length K of numeric length n_wl
#     rb_tau_k = rb_tau_k,                             # list length K of numeric length q
#     rb_L_k = rb_L_k,
#     rb_sigma2_k = as.numeric(rb_sigma2_k),           # numeric length K
#     rb_pi = as.numeric(rb_pi),                       # numeric length K (Stan: simplex[K])
#
#     # geometry
#     water_type = as.integer(water_type),
#     theta_view = as.numeric(theta_view),
#     shallow = as.integer(shallow),
#
#     # priors / hyperparameters
#     a_nap_star_mu = as.numeric(dots$a_nap_star_mu),
#     a_nap_star_sd = as.numeric(dots$a_nap_star_sd),
#     bb_p_star_mu  = as.numeric(dots$bb_p_star_mu),
#     bb_p_star_sd  = as.numeric(dots$bb_p_star_sd),
#     a_g_s_mu      = as.numeric(dots$a_g_s_mu),
#     a_g_s_sd      = as.numeric(dots$a_g_s_sd),
#     a_nap_s_mu    = as.numeric(dots$a_nap_s_mu),
#     a_nap_s_sd    = as.numeric(dots$a_nap_s_sd),
#     bb_p_gamma_mu = as.numeric(dots$bb_p_gamma_mu),
#     bb_p_gamma_sd = as.numeric(dots$bb_p_gamma_sd),
#     h_w_mu_base = as.numeric(dots$h_w_mu),
#     h_w_sd_base = as.numeric(dots$h_w_sd)
#   )
# }

#' Make base data needed by STAN model (low-rank PCA mixture bottom prior in ETA/logit space)
#'
#' @export
make_stan_data_base_gmm <- function(
    wavelength,
    water_type, theta_view, shallow,
    bottom_class_names,
    rb_obs_path = "/home/raphael/R/SABER/inst/extdata/r_b_gamache_obs.csv",
    q = 5,
    use_empirical_pi = FALSE,
    ...
) {
  dots <- list(...)

  luts <- make_luts_gmm(
    wavelength = wavelength,
    rb_obs_path = rb_obs_path,
    q = q
  )

  if (length(bottom_class_names) < 1) stop("bottom_class_names must have length >= 1")
  bad <- setdiff(bottom_class_names, luts$classes)
  if (length(bad) > 0) stop("Bottom class not found in library: ", paste(bad, collapse = ", "))
  # subset in requested order (ETA/logit-space objects)
  eta_mu_k <- lapply(bottom_class_names, function(cl) luts$eta_mu_k[[cl]])  # each length n_wl
  eta_L_k  <- lapply(bottom_class_names, function(cl) luts$eta_L_k[[cl]])   # each q x q (lower)

  sigma2_full <- setNames(luts$eta_sigma2_k, luts$classes)
  eta_sigma2_k <- unname(sigma2_full[bottom_class_names])

  # weights (empirical proportions from luts, or uniform)
  K <- length(bottom_class_names)
  if (use_empirical_pi) {
    pi_full <- setNames(luts$eta_pi, luts$classes)
    eta_pi <- unname(pi_full[bottom_class_names])
    eta_pi <- eta_pi / sum(eta_pi)
  } else {
    eta_pi <- rep(1 / K, K)
  }

  wl <- as.numeric(wavelength)
  n_wl <- length(wl)

  # hyperparameters: fail fast if missing
  hyp_needed <- c("a_nap_star_mu","a_nap_star_sd","bb_p_star_mu","bb_p_star_sd",
                  "a_g_s_mu","a_g_s_sd","a_nap_s_mu","a_nap_s_sd",
                  "bb_p_gamma_mu","bb_p_gamma_sd","h_w_mu","h_w_sd")
  hyp_miss <- setdiff(hyp_needed, names(dots))
  if (length(hyp_miss) > 0) stop("Missing required hyperparameters in `...`: ", paste(hyp_miss, collapse = ", "))

  list(
    n_wl = n_wl,
    wavelength = wl,

    a_w  = luts$a_w,
    a0   = luts$a0,
    a1   = luts$a1,
    bb_w = luts$bb_w,

    # low-rank mixture prior inputs (ETA/logit space; MUST match Stan names)
    K = as.integer(K),
    q = as.integer(luts$q),
    eta_U = luts$eta_U,                         # matrix [n_wl, q]
    eta_mu_k = eta_mu_k,                       # list length K of numeric length n_wl
    eta_L_k = eta_L_k,                         # list length K of matrix [q, q] lower
    eta_sigma2_k = as.numeric(eta_sigma2_k),   # numeric length K
    eta_pi = as.numeric(eta_pi),               # numeric length K (Stan: simplex[K])

    # geometry
    water_type = as.integer(water_type),
    theta_view = as.numeric(theta_view),
    shallow = as.integer(shallow),

    # priors / hyperparameters
    a_nap_star_mu = as.numeric(dots$a_nap_star_mu),
    a_nap_star_sd = as.numeric(dots$a_nap_star_sd),
    bb_p_star_mu  = as.numeric(dots$bb_p_star_mu),
    bb_p_star_sd  = as.numeric(dots$bb_p_star_sd),
    a_g_s_mu      = as.numeric(dots$a_g_s_mu),
    a_g_s_sd      = as.numeric(dots$a_g_s_sd),
    a_nap_s_mu    = as.numeric(dots$a_nap_s_mu),
    a_nap_s_sd    = as.numeric(dots$a_nap_s_sd),
    bb_p_gamma_mu = as.numeric(dots$bb_p_gamma_mu),
    bb_p_gamma_sd = as.numeric(dots$bb_p_gamma_sd),
    h_w_mu_base        = as.numeric(dots$h_w_mu),
    h_w_sd_base        = as.numeric(dots$h_w_sd)
  )
}

#' Make base data needed by STAN model (spline mixture bottom prior)
#'
#' @export
make_stan_data_base_splines <- function(
    wavelength,
    water_type, theta_view, shallow,
    bottom_class_names,
    rb_obs_path = "/home/raphael/R/SABER/inst/extdata/r_b_gamache_obs.csv",
    m_spline = 15,
    spline_degree = 3,
    use_empirical_pi = TRUE,
    ...
) {
  dots <- list(...)

  luts <- make_luts_gmm(
    wavelength = wavelength,
    rb_obs_path = rb_obs_path,
    q = 5, # legacy PCA still computed; irrelevant for spline model
    m_spline = m_spline,
    spline_degree = spline_degree
  )

  needed <- c("classes","K","rb_B","m_spline","rb_c_mu_k","rb_c_L_k","rb_pi",
              "a_w","a0","a1","bb_w")
  miss <- setdiff(needed, names(luts))
  if (length(miss) > 0) stop("make_luts_gmm() missing spline fields: ", paste(miss, collapse = ", "))

  if (length(bottom_class_names) < 1) stop("bottom_class_names must have length >= 1")
  bad <- setdiff(bottom_class_names, luts$classes)
  if (length(bad) > 0) stop("Bottom class not found in library: ", paste(bad, collapse = ", "))

  # subset in requested order
  rb_c_mu_k <- lapply(bottom_class_names, function(cl) luts$rb_c_mu_k[[cl]])
  rb_c_L_k  <- lapply(bottom_class_names, function(cl) luts$rb_c_L_k[[cl]])

  K <- length(bottom_class_names)
  if (use_empirical_pi) {
    pi_full <- setNames(luts$rb_pi, luts$classes)
    rb_pi <- unname(pi_full[bottom_class_names])
    rb_pi <- rb_pi / sum(rb_pi)
  } else {
    rb_pi <- rep(1 / K, K)
  }

  wl <- as.numeric(wavelength)
  n_wl <- length(wl)

  # hyperparameters: fail fast if missing
  hyp_needed <- c("a_nap_star_mu","a_nap_star_sd","bb_p_star_mu","bb_p_star_sd",
                  "a_g_s_mu","a_g_s_sd","a_nap_s_mu","a_nap_s_sd",
                  "bb_p_gamma_mu","bb_p_gamma_sd","h_w_mu","h_w_sd")
  hyp_miss <- setdiff(hyp_needed, names(dots))
  if (length(hyp_miss) > 0) stop("Missing required hyperparameters in `...`: ", paste(hyp_miss, collapse = ", "))

  list(
    n_wl = n_wl,
    wavelength = wl,

    a_w  = luts$a_w,
    a0   = luts$a0,
    a1   = luts$a1,
    bb_w = luts$bb_w,

    # spline bottom prior inputs
    K = as.integer(K),
    m = as.integer(luts$m_spline),
    rb_B = luts$rb_B,
    rb_c_mu_k = rb_c_mu_k,
    rb_c_L_k  = rb_c_L_k,
    rb_pi = as.numeric(rb_pi),

    # geometry
    water_type = as.integer(water_type),
    theta_view = as.numeric(theta_view),
    shallow = as.integer(shallow),

    # priors / hyperparameters
    a_nap_star_mu = as.numeric(dots$a_nap_star_mu),
    a_nap_star_sd = as.numeric(dots$a_nap_star_sd),
    bb_p_star_mu  = as.numeric(dots$bb_p_star_mu),
    bb_p_star_sd  = as.numeric(dots$bb_p_star_sd),
    a_g_s_mu      = as.numeric(dots$a_g_s_mu),
    a_g_s_sd      = as.numeric(dots$a_g_s_sd),
    a_nap_s_mu    = as.numeric(dots$a_nap_s_mu),
    a_nap_s_sd    = as.numeric(dots$a_nap_s_sd),
    bb_p_gamma_mu = as.numeric(dots$bb_p_gamma_mu),
    bb_p_gamma_sd = as.numeric(dots$bb_p_gamma_sd),
    h_w_mu_base = as.numeric(dots$h_w_mu),
    h_w_sd_base = as.numeric(dots$h_w_sd)
  )
}

#' Make per observation data needed by STAN model
#'
#' @export
make_stan_data_obs <- function(
    df,
    stan_data_base,
    use_measured_sigma = TRUE,
    use_depth_prior_from_obs = TRUE,
    h_w_sd_obs = NA_real_,          # allow override
    sigma_fallback = 0.00005,
    sigma_floor = 1e-8
) {
  wl <- as.numeric(df$wavelength)
  rrs_obs <- as.numeric(df$rrs_0m)
  theta_sun <- abs(unique(as.numeric(df$theta_sun)))

  stopifnot(length(wl) == stan_data_base$n_wl)
  stopifnot(all(abs(wl - stan_data_base$wavelength) < 1e-8))

  sigma_vec <- if (use_measured_sigma) as.numeric(df$rrs_unc) else rep(sigma_fallback, stan_data_base$n_wl)
  sigma_vec[!is.finite(sigma_vec)] <- NA_real_
  if (anyNA(sigma_vec)) {
    s0 <- median(sigma_vec, na.rm = TRUE)
    sigma_vec[is.na(sigma_vec)] <- if (is.finite(s0)) s0 else sigma_fallback
  }
  sigma_vec <- pmax(sigma_vec, sigma_floor)

  # depth info from df (if present)
  h_w <- abs(unique(as.numeric(df$h_w)))
  if (length(h_w) != 1 || !is.finite(h_w)) h_w <- NA_real_

  out <- list(
    rrs_obs = rrs_obs,
    sigma_rrs = sigma_vec,
    theta_sun = theta_sun,
    wl = wl
  )

  if (use_depth_prior_from_obs) {
    out$h_w_mu_obs <- h_w
    out$h_w_sd_obs <- as.numeric(h_w_sd_obs)
  }

  out
}

#' Make stan data
#'
#' @export
make_stan_data <- function(
    stan_data_base,
    stan_data_obs,
    h_w_mu_from_obs = F,
    h_w_sd_from_obs = F
    ) {

  if (h_w_mu_from_obs) {
    h_w_mu <- stan_data_obs$h_w_mu_obs
  } else {
    h_w_mu <- stan_data_base$h_w_mu_base
  }

  if (h_w_sd_from_obs) {
    h_w_sd <- stan_data_obs$h_w_sd_obs
  } else {
    h_w_sd <- stan_data_base$h_w_sd_base
  }

  out <- append(stan_data_base, stan_data_obs)
  out$h_w_mu_obs <- NULL
  out$h_w_mu_base <- NULL
  out$h_w_sd_obs <- NULL
  out$h_w_sd_base <- NULL
  out$h_w_mu <- h_w_mu
  out$h_w_sd <- h_w_sd

  return(out)
}
