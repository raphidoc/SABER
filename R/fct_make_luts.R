interp1 <- function(x, y, xout) {
  approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
}

#' @import dplyr
#'
#' @export
make_luts <- function(
    wavelength,
    sd_floor = 1e-6
) {

  stopifnot(is.numeric(wavelength), length(wavelength) >= 2)

  # data("a_w", package = "SABER", envir = environment())
  # data("a0_a1_phyto", package = "SABER", envir = environment())
  # data("r_b_gamache", package = "SABER", envir = environment())

  a_w_int  <- interp1(a_w$wavelength,  a_w$a_w, wavelength)
  a0_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a0, wavelength)
  a1_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a1, wavelength)

  compute_bb_w <- function(wl) 0.00111 * (wl / 500)^(-4.32)
  bb_w <- compute_bb_w(wavelength)

  rb <- r_b_gamache %>%
    transmute(
      class = as.character(.data$class),
      wl    = .data$wavelength,
      mu    = .data$r_b_mean,
      sd    = .data$r_b_sd
    )

  classes <- rb %>% distinct(class) %>% arrange(class) %>% pull(class)

  # Wide mean + sd tables on the library's native wavelength grid
  rb_mu_wide <- rb %>%
    select(class, wl, mu) %>%
    pivot_wider(names_from = class, values_from = mu) %>%
    arrange(wl)

  rb_sd_wide <- rb %>%
    select(class, wl, sd) %>%
    pivot_wider(names_from = class, values_from = sd) %>%
    arrange(wl)

  rb_wl_master <- rb_mu_wide$wl

  # Interpolate each class mean and sd onto the Stan wavelength grid
  r_b_mu_lib <- vapply(classes, function(cl) {
    interp1(rb_wl_master, rb_mu_wide[[cl]], wavelength)
  }, FUN.VALUE = numeric(length(wavelength)))

  r_b_sd_lib <- vapply(classes, function(cl) {
    interp1(rb_wl_master, rb_sd_wide[[cl]], wavelength)
  }, FUN.VALUE = numeric(length(wavelength)))

  r_b_mu_lib <- matrix(r_b_mu_lib, nrow = length(wavelength), ncol = length(classes))
  r_b_sd_lib <- matrix(r_b_sd_lib, nrow = length(wavelength), ncol = length(classes))

  colnames(r_b_mu_lib) <- classes
  colnames(r_b_sd_lib) <- classes

  # Safety: enforce nonnegative sd, add floor to avoid zeros in Stan
  r_b_sd_lib <- pmax(r_b_sd_lib, sd_floor)

  list(
    a_w = a_w_int, a0 = a0_int, a1 = a1_int, bb_w = bb_w,
    r_b_mu_lib = r_b_mu_lib,
    r_b_sd_lib = r_b_sd_lib,
    classes = classes
  )
}

# #' @import dplyr
# #' @export
# make_luts_gmm <- function(
#     wavelength,
#     sd_floor = 1e-6,
#     rb_obs_path = "/home/raphael/R/SABER/inst/extdata/r_b_gamache_obs.csv",
#     q = 5,
#     sigma2_floor = 1e-6,
#     tau_floor = 1e-6#,
#     # m_spline = 15,
#     # spline_degree = 3,
#     # eps = 1e-6,
#     # ridge = 1e-8,
#     # cov_jitter = 1e-6
# ) {
#   stopifnot(is.numeric(wavelength), length(wavelength) >= 2)
#
#   # ---- IOP / water LUTs ----
#   a_w_int  <- interp1(a_w$wavelength,  a_w$a_w, wavelength)
#   a0_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a0, wavelength)
#   a1_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a1, wavelength)
#
#   compute_bb_w <- function(wl) 0.00111 * (wl / 500)^(-4.32)
#   bb_w <- compute_bb_w(wavelength)
#
#   # ---- Bottom library: mean/sd table (legacy; keep) ----
#   rb_ms <- r_b_gamache %>%
#     dplyr::transmute(
#       class = as.character(.data$class),
#       wl    = as.numeric(.data$wavelength),
#       mu    = as.numeric(.data$r_b_mean),
#       sd    = as.numeric(.data$r_b_sd)
#     )
#
#   classes_ms <- rb_ms %>%
#     dplyr::distinct(class) %>%
#     dplyr::arrange(class) %>%
#     dplyr::pull(class)
#
#   rb_mu_wide <- rb_ms %>%
#     dplyr::select(class, wl, mu) %>%
#     tidyr::pivot_wider(names_from = class, values_from = mu) %>%
#     dplyr::arrange(wl)
#
#   rb_sd_wide <- rb_ms %>%
#     dplyr::select(class, wl, sd) %>%
#     tidyr::pivot_wider(names_from = class, values_from = sd) %>%
#     dplyr::arrange(wl)
#
#   rb_wl_master <- rb_mu_wide$wl
#
#   r_b_mu_lib <- vapply(classes_ms, function(cl) {
#     interp1(rb_wl_master, rb_mu_wide[[cl]], wavelength)
#   }, FUN.VALUE = numeric(length(wavelength)))
#
#   r_b_sd_lib <- vapply(classes_ms, function(cl) {
#     interp1(rb_wl_master, rb_sd_wide[[cl]], wavelength)
#   }, FUN.VALUE = numeric(length(wavelength)))
#
#   r_b_mu_lib <- matrix(r_b_mu_lib, nrow = length(wavelength), ncol = length(classes_ms))
#   r_b_sd_lib <- matrix(r_b_sd_lib, nrow = length(wavelength), ncol = length(classes_ms))
#   colnames(r_b_mu_lib) <- classes_ms
#   colnames(r_b_sd_lib) <- classes_ms
#   r_b_sd_lib <- pmax(r_b_sd_lib, sd_floor)
#
#   # ---- Bottom replicate spectra: low-rank PCA prior ----
#   rb_obs <- readr::read_csv(rb_obs_path, show_col_types = FALSE) %>%
#     dplyr::mutate(
#       class = as.character(.data$class),
#       wavelength = as.numeric(.data$wavelength),
#       r_b = as.numeric(.data$r_b)
#     )
#
#   # identify replicate id if present; otherwise fabricate within class
#   id_col <- intersect(names(rb_obs), c("time_target"))
#   if (length(id_col) == 0) {
#     # rb_obs <- rb_obs %>%
#     #   group_by(class) %>%
#     #   arrange(wavelength, .by_group = TRUE) %>%
#     #   mutate(rep_id = cumsum(c(1, diff(wavelength) < 0))) %>%
#     #   ungroup()
#     # id_col <- "rep_id"
#   } else {
#     id_col <- id_col[1]
#   }
#
#   spectra_df <- rb_obs %>%
#     dplyr::select(class, dplyr::all_of(id_col), wavelength, r_b) %>%
#     tidyr::pivot_wider(names_from = wavelength, values_from = r_b, names_prefix = "wl_")
#
#   wl_cols <- grep("^wl_", names(spectra_df), value = TRUE)
#   wl_native <- as.numeric(sub("^wl_", "", wl_cols))
#   ord <- order(wl_native)
#   wl_cols <- wl_cols[ord]
#   wl_native <- wl_native[ord]
#
#   X_native <- as.matrix(spectra_df[, wl_cols, drop = FALSE])
#   wl_target <- as.numeric(wavelength)
#
#   # interpolate each replicate spectrum onto Stan grid
#   X <- t(apply(X_native, 1, function(row_vals) {
#     y <- as.numeric(row_vals)
#     ok <- is.finite(y)
#     if (sum(ok) < 2) return(rep(NA_real_, length(wl_target)))
#     interp1(wl_native[ok], y[ok], wl_target)
#   }))  # [n_spec, n_wl]
#
#   # fill any non-finite with column medians
#   if (any(!is.finite(X))) {
#     col_med <- apply(X, 2, function(z) stats::median(z, na.rm = TRUE))
#     for (j in seq_len(ncol(X))) {
#       bad <- !is.finite(X[, j])
#       if (any(bad)) X[bad, j] <- col_med[j]
#     }
#   }
#
#   cls_vec <- spectra_df$class
#   classes <- sort(unique(cls_vec))
#   K <- length(classes)
#
#   eps <- 1e-6
#   clamp01 <- function(x) pmin(pmax(x, eps), 1 - eps)
#
#   eta_X <- qlogis(clamp01(X))  # X is your [n_spec, n_wl] reflectance matrix
#
#   mu_global <- colMeans(eta_X)
#   Xc <- sweep(eta_X, 2, mu_global, "-")
#   p <- prcomp(Xc, center = FALSE, scale. = FALSE)
#   q_eff <- min(q, ncol(p$rotation), ncol(X))
#   eta_U <- p$rotation[, 1:q_eff, drop=FALSE]
#
#   # global PCA basis (shared across classes)
#   # mu_global <- colMeans(X)
#   # Xc <- sweep(X, 2, mu_global, "-")
#   # p <- stats::prcomp(Xc, center = FALSE, scale. = FALSE)
#   # q_eff <- min(q, ncol(p$rotation), ncol(X))
#   # rb_U <- p$rotation[, 1:q_eff, drop = FALSE]  # [n_wl, q]
#
#   rb_mu_k <- vector("list", K)
#   rb_tau_k <- vector("list", K)
#   rb_L_k <- vector("list", K)
#   rb_sigma2_k <- numeric(K)
#   counts <- integer(K)
#
#   for (i in seq_along(classes)) {
#     cl <- classes[i]
#     idx <- which(cls_vec == cl)
#     counts[i] <- length(idx)
#
#     Xk <- X[idx, , drop = FALSE]
#     mu_k <- colMeans(Xk)
#     rb_mu_k[[i]] <- as.numeric(mu_k)
#
#     # residuals around class mean
#     # Rk <- sweep(Xk, 2, mu_k, "-")        # [n_k, n_wl]
#     # Ak <- Rk %*% rb_U                    # [n_k, q]
#     # tau <- apply(Ak, 2, stats::sd)
#     # tau[!is.finite(tau)] <- 0
#     # tau <- pmax(tau, tau_floor)
#     # rb_tau_k[[i]] <- as.numeric(tau)
#     #
#     # Sk <- stats::cov(Ak)                         # q x q
#     # # jitter for numerical stability
#     # Sk <- Sk + diag(tau_floor^2, nrow(Sk))
#     # Lk <- chol(Sk)                               # upper by default in R
#     # rb_L_k[[i]] <- t(Lk)                         # convert to lower for Stan
#     #
#     # Rhat <- Ak %*% t(rb_U)
#     # Ek <- Rk - Rhat
#     # sigma2 <- mean(Ek^2)
#     # sigma2 <- max(sigma2, sigma2_floor)
#     # rb_sigma2_k[i] <- sigma2
#
#     eta_mu_k[[i]] <- colMeans(eta_X[idx, , drop=FALSE])
#
#     Rk <- sweep(eta_X[idx, , drop=FALSE], 2, eta_mu_k[[i]], "-")
#     Ak <- Rk %*% eta_U
#     Sk <- cov(Ak) + diag(jitter, ncol(Ak))
#     eta_L_k[[i]] <- t(chol(Sk))     # lower for Stan
#     eta_sigma2_k[i] <- max(mean((Rk - Ak %*% t(eta_U))^2), sigma2_floor)
#
#   }
#
#   names(rb_mu_k) <- classes
#   names(rb_tau_k) <- classes
#   names(rb_L_k) <- classes
#   rb_pi <- counts / sum(counts)
#
#   # # ---- Bottom replicate spectra: spline prior in logit space (NEW) ----
#   # # Build spline basis on the Stan wavelength grid. Include intercept.
#   # B <- splines::bs(
#   #   wl_target,
#   #   df = m_spline,
#   #   degree = spline_degree,
#   #   intercept = TRUE
#   # )
#   # B <- as.matrix(B)
#   # m_eff <- ncol(B)
#   #
#   # # Fit spline coefficients for each replicate in logit space:
#   # # eta(λ) = B %*% c,  r_b(λ) = inv_logit(eta)
#   # clamp01 <- function(x, eps) pmin(pmax(x, eps), 1 - eps)
#   #
#   # BtB <- crossprod(B)  # m x m
#   # BtB_ridge <- BtB + diag(ridge, m_eff)
#   # BtB_inv <- solve(BtB_ridge)
#   #
#   # coef_mat <- matrix(NA_real_, nrow = nrow(X), ncol = m_eff)
#   # for (i in seq_len(nrow(X))) {
#   #   rb_i <- clamp01(as.numeric(X[i, ]), eps)
#   #   eta_i <- qlogis(rb_i)
#   #   coef_mat[i, ] <- BtB_inv %*% crossprod(B, eta_i)
#   # }
#   #
#   # # Class-conditional mean + covariance in coefficient space
#   # rb_c_mu_k <- vector("list", K)
#   # rb_c_L_k  <- vector("list", K)
#   #
#   # for (i in seq_along(classes)) {
#   #   cl <- classes[i]
#   #   idx <- which(cls_vec == cl)
#   #
#   #   Ck <- coef_mat[idx, , drop = FALSE]  # n_k x m
#   #   mu_c <- colMeans(Ck)
#   #   rb_c_mu_k[[i]] <- as.numeric(mu_c)
#   #
#   #   Sk <- stats::cov(Ck)
#   #   if (!all(is.finite(Sk))) {
#   #     Sk <- diag(1, m_eff) * cov_jitter
#   #   }
#   #   Sk <- Sk + diag(cov_jitter, m_eff)
#   #
#   #   # R chol is upper; Stan wants lower
#   #   Lk <- chol(Sk)
#   #   rb_c_L_k[[i]] <- t(Lk)
#   # }
#   #
#   # names(rb_c_mu_k) <- classes
#   # names(rb_c_L_k)  <- classes
#
#   list(
#     # IOP LUTs
#     a_w = a_w_int, a0 = a0_int, a1 = a1_int, bb_w = bb_w,
#
#     # legacy mean/sd (optional)
#     r_b_mu_lib = r_b_mu_lib,
#     r_b_sd_lib = r_b_sd_lib,
#
#     # low-rank PCA prior objects
#     classes = classes,
#     K = as.integer(K),
#     q = as.integer(q_eff),
#     rb_U = rb_U,
#     rb_mu_k = rb_mu_k,
#     rb_tau_k = rb_tau_k,
#     rb_L_k = rb_L_k,
#     rb_sigma2_k = as.numeric(rb_sigma2_k),
#     rb_pi = as.numeric(rb_pi)#,
#     # ---- spline bottom prior objects (NEW) ----
#     # rb_B = B,                       # [n_wl, m]
#     # m_spline = as.integer(m_eff),
#     # rb_c_mu_k = rb_c_mu_k,          # list length K, each length m
#     # rb_c_L_k  = rb_c_L_k            # list length K, each m x m (lower)
#   )
# }

#' @import dplyr
#' @export
make_luts_gmm <- function(
    wavelength,
    sd_floor = 1e-6,
    rb_obs_path = "/home/raphael/R/SABER/inst/extdata/r_b_gamache_obs.csv",
    q = 5,
    sigma2_floor = 1e-6,
    eps = 1e-6,
    cov_jitter = 1e-6
) {
  stopifnot(is.numeric(wavelength), length(wavelength) >= 2)

  # ---- IOP / water LUTs ----
  a_w_int  <- interp1(a_w$wavelength,  a_w$a_w, wavelength)
  a0_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a0, wavelength)
  a1_int   <- interp1(a0_a1_phyto$wavelength, a0_a1_phyto$a1, wavelength)

  compute_bb_w <- function(wl) 0.00111 * (wl / 500)^(-4.32)
  bb_w <- compute_bb_w(wavelength)

  # ---- Bottom library: mean/sd table (legacy; keep) ----
  rb_ms <- r_b_gamache %>%
    dplyr::transmute(
      class = as.character(.data$class),
      wl    = as.numeric(.data$wavelength),
      mu    = as.numeric(.data$r_b_mean),
      sd    = as.numeric(.data$r_b_sd)
    )

  classes_ms <- rb_ms %>%
    dplyr::distinct(class) %>%
    dplyr::arrange(class) %>%
    dplyr::pull(class)

  rb_mu_wide <- rb_ms %>%
    dplyr::select(class, wl, mu) %>%
    tidyr::pivot_wider(names_from = class, values_from = mu) %>%
    dplyr::arrange(wl)

  rb_sd_wide <- rb_ms %>%
    dplyr::select(class, wl, sd) %>%
    tidyr::pivot_wider(names_from = class, values_from = sd) %>%
    dplyr::arrange(wl)

  rb_wl_master <- rb_mu_wide$wl

  r_b_mu_lib <- vapply(classes_ms, function(cl) {
    interp1(rb_wl_master, rb_mu_wide[[cl]], wavelength)
  }, FUN.VALUE = numeric(length(wavelength)))

  r_b_sd_lib <- vapply(classes_ms, function(cl) {
    interp1(rb_wl_master, rb_sd_wide[[cl]], wavelength)
  }, FUN.VALUE = numeric(length(wavelength)))

  r_b_mu_lib <- matrix(r_b_mu_lib, nrow = length(wavelength), ncol = length(classes_ms))
  r_b_sd_lib <- matrix(r_b_sd_lib, nrow = length(wavelength), ncol = length(classes_ms))
  colnames(r_b_mu_lib) <- classes_ms
  colnames(r_b_sd_lib) <- classes_ms
  r_b_sd_lib <- pmax(r_b_sd_lib, sd_floor)

  # ---- Bottom replicate spectra: LOW-RANK PCA PRIOR IN LOGIT SPACE ----
  rb_obs <- readr::read_csv(rb_obs_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      class = as.character(.data$class),
      wavelength = as.numeric(.data$wavelength),
      r_b = as.numeric(.data$r_b)
    )

  # replicate identifier (required)
  id_col <- intersect(names(rb_obs), c("time_target"))
  if (length(id_col) == 0) {
    stop("No replicate identifier column found (expected 'time_target'). Add a replicate id column or adapt this function.")
  }
  id_col <- id_col[1]

  spectra_df <- rb_obs %>%
    dplyr::select(class, dplyr::all_of(id_col), wavelength, r_b) %>%
    tidyr::pivot_wider(names_from = wavelength, values_from = r_b, names_prefix = "wl_")

  wl_cols <- grep("^wl_", names(spectra_df), value = TRUE)
  wl_native <- as.numeric(sub("^wl_", "", wl_cols))
  ord <- order(wl_native)
  wl_cols <- wl_cols[ord]
  wl_native <- wl_native[ord]

  X_native <- as.matrix(spectra_df[, wl_cols, drop = FALSE])
  wl_target <- as.numeric(wavelength)

  # interpolate each replicate spectrum onto Stan grid
  X <- t(apply(X_native, 1, function(row_vals) {
    y <- as.numeric(row_vals)
    ok <- is.finite(y)
    if (sum(ok) < 2) return(rep(NA_real_, length(wl_target)))
    interp1(wl_native[ok], y[ok], wl_target)
  }))  # [n_spec, n_wl]

  # fill any non-finite with column medians
  if (any(!is.finite(X))) {
    col_med <- apply(X, 2, function(z) stats::median(z, na.rm = TRUE))
    for (j in seq_len(ncol(X))) {
      bad <- !is.finite(X[, j])
      if (any(bad)) X[bad, j] <- col_med[j]
    }
  }

  cls_vec <- spectra_df$class
  classes <- sort(unique(cls_vec))
  K <- length(classes)

  clamp01 <- function(x) pmin(pmax(x, eps), 1 - eps)

  # reflectance -> unconstrained latent eta (logit space)
  eta_X <- qlogis(clamp01(X))  # [n_spec, n_wl]

  # global PCA basis (shared across classes) in eta-space
  mu_global <- colMeans(eta_X)
  Xc <- sweep(eta_X, 2, mu_global, "-")
  p <- stats::prcomp(Xc, center = FALSE, scale. = FALSE)

  q_eff <- min(q, ncol(p$rotation), ncol(eta_X))
  eta_U <- p$rotation[, 1:q_eff, drop = FALSE]  # [n_wl, q]

  # class-conditional stats in eta-space
  eta_mu_k <- vector("list", K)
  eta_L_k <- vector("list", K)
  eta_sigma2_k <- numeric(K)
  counts <- integer(K)

  for (i in seq_along(classes)) {
    cl <- classes[i]
    idx <- which(cls_vec == cl)
    counts[i] <- length(idx)

    mu_k <- colMeans(eta_X[idx, , drop = FALSE])
    eta_mu_k[[i]] <- as.numeric(mu_k)

    Rk <- sweep(eta_X[idx, , drop = FALSE], 2, mu_k, "-")  # [n_k, n_wl]
    Ak <- Rk %*% eta_U                                     # [n_k, q]

    Sk <- stats::cov(Ak)
    if (!all(is.finite(Sk))) Sk <- diag(1, q_eff) * cov_jitter
    Sk <- Sk + diag(cov_jitter, q_eff)

    # R chol is upper; Stan wants lower
    eta_L_k[[i]] <- t(chol(Sk))

    Rhat <- Ak %*% t(eta_U)
    Ek <- Rk - Rhat
    eta_sigma2_k[i] <- max(mean(Ek^2), sigma2_floor)
  }

  names(eta_mu_k) <- classes
  names(eta_L_k) <- classes
  eta_pi <- counts / sum(counts)

  list(
    # IOP LUTs
    a_w = a_w_int, a0 = a0_int, a1 = a1_int, bb_w = bb_w,

    # legacy mean/sd (optional)
    r_b_mu_lib = r_b_mu_lib,
    r_b_sd_lib = r_b_sd_lib,

    # low-rank PCA prior objects (ETA / LOGIT SPACE)
    classes = classes,
    K = as.integer(K),
    q = as.integer(q_eff),

    # keep Stan names if you want minimal Stan edits
    eta_U = eta_U,
    eta_mu_k = eta_mu_k,
    eta_L_k = eta_L_k,
    eta_sigma2_k = as.numeric(eta_sigma2_k),
    eta_pi = as.numeric(eta_pi),

    # helpful metadata
    eps = eps,
    cov_jitter = cov_jitter
  )
}
