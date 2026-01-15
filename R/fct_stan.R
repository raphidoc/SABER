

# data builders -----------------------------------------------------------

interp1 <- function(x, y, xout) {
  approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
}

build_saber_stan_luts <- function(wavelength,
                                  a_w_df,
                                  a0_a1_df,
                                  r_b_long,
                                  bb_w_fun = function(wl) 0.00111 * (wl / 500)^(-4.32),
                                  sd_floor = 1e-6) {

  stopifnot(is.numeric(wavelength), length(wavelength) >= 2)

  a_w  <- interp1(a_w_df$wavelength,  a_w_df$a_w, wavelength)
  a0   <- interp1(a0_a1_df$wavelength, a0_a1_df$a0, wavelength)
  a1   <- interp1(a0_a1_df$wavelength, a0_a1_df$a1, wavelength)
  bb_w <- bb_w_fun(wavelength)

  rb <- r_b_long %>%
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

  # Spectral damping for SD: 1.0 from 350–700 nm, 0.3 above 700 nm
  sd_coef <- ifelse(wavelength <= 700, 1.0, 0.25)
  r_b_sd_lib <- r_b_sd_lib * sd_coef

  # Safety: enforce nonnegative sd, add floor to avoid zeros in Stan
  r_b_sd_lib <- pmax(r_b_sd_lib, sd_floor)

  list(
    a_w = a_w, a0 = a0, a1 = a1, bb_w = bb_w,
    r_b_mu_lib = r_b_mu_lib,
    r_b_sd_lib = r_b_sd_lib,
    classes = classes
  )
}

# make_saber_stan_data_base <- function(wavelength,
#                                       water_type, theta_sun_deg, theta_view_deg, shallow,
#                                       bottom_class_names,
#                                       pkgname = "SABER") {
#
#   data("a_w", package = pkgname, envir = environment())
#   data("a0_a1_phyto", package = pkgname, envir = environment())
#   data("r_b_gamache", package = pkgname, envir = environment())
#
#   luts <- build_saber_stan_luts(
#     wavelength = wavelength,
#     a_w_df = a_w,
#     a0_a1_df = a0_a1_phyto,
#     r_b_long = r_b_gamache
#   )
#
#   bottom_class_ids <- match(bottom_class_names, luts$classes)
#   if (any(is.na(bottom_class_ids))) {
#     missing <- bottom_class_names[is.na(bottom_class_ids)]
#     stop("Bottom class not found in library: ", paste(missing, collapse = ", "))
#   }
#
#   list(
#     n_wl = length(wavelength),
#     wavelength = as.numeric(wavelength),
#
#     # placeholders to fill per spectrum
#     rrs_obs   = rep(NA_real_, length(wavelength)),
#     sigma_rrs = rep(NA_real_, length(wavelength)),
#
#     a_w = luts$a_w,
#     a0  = luts$a0,
#     a1  = luts$a1,
#     bb_w = luts$bb_w,
#
#     n_class = length(luts$classes),
#
#     # NEW names expected by the adapted Stan model
#     r_b_mu_lib = luts$r_b_mu_lib,
#     r_b_sd_lib = luts$r_b_sd_lib,
#
#     bottom_class_ids = as.integer(bottom_class_ids),
#
#     water_type = as.integer(water_type),
#     theta_sun_deg = as.numeric(theta_sun_deg),
#     theta_view_deg = as.numeric(theta_view_deg),
#     shallow = as.integer(shallow)
#   )
# }


make_saber_stan_data_base <- function(wavelength,
                                      water_type, theta_sun_deg, theta_view_deg, shallow,
                                      bottom_class_names,
                                      K = 3,
                                      delta_scale = 0.5,
                                      basis = c("cosine", "poly"),
                                      pkgname = "SABER") {
  basis <- match.arg(basis)

  data("a_w", package = pkgname, envir = environment())
  data("a0_a1_phyto", package = pkgname, envir = environment())
  data("r_b_gamache", package = pkgname, envir = environment())

  luts <- build_saber_stan_luts(
    wavelength = wavelength,
    a_w_df = a_w,
    a0_a1_df = a0_a1_phyto,
    r_b_long = r_b_gamache
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

  if (!is.numeric(K) || length(K) != 1 || K < 1) stop("K must be a single integer >= 1.")
  K <- as.integer(K)

  if (!is.numeric(delta_scale) || length(delta_scale) != 1 || delta_scale <= 0) {
    stop("delta_scale must be a single numeric > 0.")
  }
  delta_scale <- as.numeric(delta_scale)

  # normalize wavelength to [0,1] for stable basis
  t <- (wl - min(wl)) / (max(wl) - min(wl))
  if (basis == "cosine") {
    # columns: cos(pi * k * t), k = 1..K (no intercept; amplitude already handles global level)
    Phi <- sapply(seq_len(K), function(k) cos(pi * k * t))
  } else {
    # polynomial basis: t^k, k = 1..K
    Phi <- sapply(seq_len(K), function(k) t^k)
  }
  Phi <- matrix(as.numeric(Phi), nrow = n_wl, ncol = K)

  list(
    n_wl = n_wl,
    wavelength = wl,

    # placeholders to fill per spectrum
    rrs_obs   = rep(NA_real_, n_wl),
    sigma_rrs = rep(NA_real_, n_wl),

    a_w  = luts$a_w,
    a0   = luts$a0,
    a1   = luts$a1,
    bb_w = luts$bb_w,

    n_class = length(luts$classes),

    r_b_mu_lib = luts$r_b_mu_lib,
    r_b_sd_lib = luts$r_b_sd_lib,

    bottom_class_ids = bottom_class_ids,

    # NEW for the SD-constrained low-rank deviation
    K = K,
    Phi = Phi,
    delta_scale = delta_scale,

    water_type = as.integer(water_type),
    theta_sun_deg = as.numeric(theta_sun_deg),
    theta_view_deg = as.numeric(theta_view_deg),
    shallow = as.integer(shallow)
  )
}

prepare_obs_inputs <- function(df, stan_data_base, use_measured_sigma = TRUE,
                               sigma_fallback = 0.00005, sigma_floor = 1e-8) {

  wl <- as.numeric(df$wavelength)
  rrs_obs <- as.numeric(df$rrs_0m)
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

  list(rrs_obs = rrs_obs, sigma_rrs = sigma_vec, h_w = h_w, wl = wl)
}

# stan_data_from_obs <- function(stan_data_base, rrs_obs, sigma_vec) {
#   stan_data_base$rrs_obs <- rrs_obs
#   stan_data_base$sigma_rrs   <- sigma_vec
#   stan_data_base
# }

stan_data_from_obs <- function(stan_data_base, rrs_obs, sigma_vec, h_w) {
  stan_data_base$rrs_obs <- rrs_obs
  stan_data_base$sigma_rrs   <- sigma_vec
  stan_data_base$h_w <- h_w
  stan_data_base
}



# sampling ----------------------------------------------------------------

fit_one <- function(mod, stan_data,
                    mode = c("sample", "optimize"),
                    seed = 123,
                    chains = 2, iter_warmup = 300, iter_sampling = 300,
                    parallel_chains = chains,
                    adapt_delta = 0.85, max_treedepth = 12,
                    threads_per_chain = 1,
                    output_dir = NULL) {

  mode <- match.arg(mode)

  if (mode == "optimize") {
    return(mod$optimize(
      data = stan_data,
      seed = seed,
      output_dir = output_dir,
      refresh = 0
    ))
  }

  mod$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    threads_per_chain = threads_per_chain,
    refresh = 0,
    output_dir = output_dir
  )
}

run_one_ensemble <- function(df, ensemble_id, mod, stan_data_base,
                             out_dir,
                             use_measured_sigma = TRUE,
                             mode = c("sample", "optimize"),
                             seed = 123,
                             chains = 2, iter_warmup = 300, iter_sampling = 300,
                             adapt_delta = 0.85, max_treedepth = 12,
                             threads_per_chain = 1,
                             spec_probs = c(0.05, 0.95),
                             keep_fit_rds = FALSE) {

  mode <- match.arg(mode)
  stopifnot(dir.exists(out_dir))

  # obs <- prepare_obs_inputs(df, stan_data_base, use_measured_sigma = use_measured_sigma)
  # stan_data <- stan_data_from_obs(stan_data_base, obs$rrs_obs, obs$sigma_rrs)

  obs <- prepare_obs_inputs(df, stan_data_base, use_measured_sigma = use_measured_sigma)
  stan_data <- stan_data_from_obs(stan_data_base, obs$rrs_obs, obs$sigma_rrs, obs$h_w)

  ens_dir <- file.path(out_dir, paste0("ensemble_", ensemble_id))
  dir.create(ens_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- fit_one(
    mod = mod,
    stan_data = stan_data,
    mode = mode,
    seed = seed,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    threads_per_chain = threads_per_chain,
    output_dir = ens_dir
  )

  # Save inputs for reproducibility (small)
  saveRDS(list(ensemble = ensemble_id,
               wavelength = obs$wl,
               sigma_rrs = obs$sigma_rrs),
          file = file.path(ens_dir, "inputs_used.rds"))

  if (mode == "optimize") {
    # MAP outputs are different; keep it light
    saveRDS(fit, file = file.path(ens_dir, "fit_map.rds"))
    return(list(ensemble = ensemble_id, out_dir = ens_dir, mode = "optimize"))
  }

  # Diagnostics
  diag_samp <- sampler_diag_stats(fit)
  summ_all  <- fit$summary()
  diag_fit  <- tibble(
    max_rhat = max(summ_all$rhat, na.rm = TRUE),
    min_ess_bulk = min(summ_all$ess_bulk, na.rm = TRUE)
  )
  diag <- bind_cols(diag_samp, diag_fit)
  saveRDS(diag, file = file.path(ens_dir, "diagnostics.rds"))

  # Params
  par_sum <- summarize_params_fast(fit)
  saveRDS(par_sum, file = file.path(ens_dir, "param_summary.rds"))

  par_mean_wide <- par_sum %>%
    # filter(
    #   variable %in% c("chl","a_g_440","a_nap_440","bb_p_550","h_w") |
    #     grepl("^r_b*", variable)
    # ) %>%
    filter(
      variable %in% c("chl","a_g_440","a_nap_440","bb_p_550") |
        grepl("^r_b*", variable)
    ) %>%
    select(variable, mean) %>%
    pivot_wider(names_from = variable, values_from = mean)
  saveRDS(par_mean_wide, file = file.path(ens_dir, "param_mean_wide.rds"))

  # Spectrum summaries + metrics (fast quantiles)
  spec <- summarize_rrs_hat_fast(fit, probs = spec_probs)
  metrics <- fit_metrics_1(obs$rrs_obs, spec$rrs_hat, spec$rrs_q_lo, spec$rrs_q_hi, sigma_vec = obs$sigma_rrs)

  rrs_out <- list(
    wl = obs$wl,
    mean = spec$rrs_hat,
    q_lo = spec$rrs_q_lo,
    q_hi = spec$rrs_q_hi,
    probs = spec_probs,
    metrics = metrics
  )
  saveRDS(rrs_out, file = file.path(ens_dir, "rrs_posterior_summaries.rds"))

  # Optional: save fit object (usually unnecessary; CSVs are already in ens_dir)
  if (keep_fit_rds) saveRDS(fit, file = file.path(ens_dir, "fit_object.rds"))

  list(
    ensemble = ensemble_id,
    out_dir = ens_dir,
    diag = diag,
    par_mean = par_mean_wide,
    par_summary = par_sum,
    rrs = rrs_out
  )
}


# post-processing ---------------------------------------------------------

summarize_rrs_hat_fast <- function(fit, probs = c(0.05, 0.95)) {
  rrs_mat <- posterior::as_draws_matrix(fit$draws("rrs_hat"))  # draws x n_wl
  rrs_hat <- colMeans(rrs_mat)
  rrs_qs <- matrixStats::colQuantiles(rrs_mat, probs = probs, drop = FALSE)

  # r_b_mat <- posterior::as_draws_matrix(fit$draws("r_b_hat"))  # draws x n_wl
  # r_b_hat <- colMeans(r_b_mat)
  # r_b_qs <- matrixStats::colQuantiles(r_b_mat, probs = probs, drop = FALSE)

  tibble(
    rrs_hat = rrs_hat,
    rrs_q_lo = rrs_qs[,1],
    rrs_q_hi = rrs_qs[,2]#,
    # r_b_hat = r_b_hat,
    # r_b_q_lo = r_b_qs[,1],
    # r_b_q_hi = r_b_qs[,2],
  )
}

fit_metrics_1 <- function(rrs_obs, mean_hat, q_lo, q_hi, sigma_vec = NULL) {
  resid <- rrs_obs - mean_hat
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))

  dotp <- sum(rrs_obs * mean_hat, na.rm = TRUE)
  n1   <- sqrt(sum(rrs_obs^2, na.rm = TRUE))
  n2   <- sqrt(sum(mean_hat^2, na.rm = TRUE))
  sam  <- acos(max(-1, min(1, dotp / (n1 * n2))))

  coverage <- mean(rrs_obs >= q_lo & rrs_obs <= q_hi, na.rm = TRUE)

  z_rmse <- NA_real_
  if (!is.null(sigma_vec)) {
    z <- (rrs_obs - mean_hat) / sigma_vec
    z_rmse <- sqrt(mean(z^2, na.rm = TRUE))
  }

  tibble(rmse = rmse, sam = sam, coverage = coverage, z_rmse = z_rmse)
}

summarize_params_fast <- function(fit,
                                  par_vars =
                                    # c("chl","a_g_440","a_nap_440","bb_p_550","h_w","bottom_mix","a_bottom"),
                                    # c("chl","a_g_440","a_nap_440","bb_p_550","h_w","bottom_mix"),
                                    # c("chl","a_g_440","a_nap_440", "a_nap_s", "a_g_s", "bb_p_550", "bb_p_gamma", "h_w","r_b_mix", "r_b_a", "sigma_model"),
                                    c("chl","a_g_440","a_nap_440", "a_nap_s", "a_g_s", "bb_p_550", "bb_p_gamma", "r_b_mix", "r_b_a", "sigma_model"),
                                  probs = c(0.05, 0.5, 0.95)) {
  posterior::summarise_draws(
    fit$draws(variables = par_vars),
    mean, sd,
    ~posterior::quantile2(.x, probs = probs),
    rhat, ess_bulk, ess_tail
  )
}

sampler_diag_stats <- function(fit) {
  d <- posterior::as_draws_matrix(fit$sampler_diagnostics())
  tibble(
    n_divergent = if ("divergent__" %in% colnames(d)) sum(d[, "divergent__"]) else NA_integer_,
    max_treedepth = if ("treedepth__" %in% colnames(d)) max(d[, "treedepth__"]) else NA_integer_,
    mean_treedepth = if ("treedepth__" %in% colnames(d)) mean(d[, "treedepth__"]) else NA_real_,
    mean_stepsize = if ("stepsize__" %in% colnames(d)) mean(d[, "stepsize__"]) else NA_real_,
    mean_n_leapfrog = if ("n_leapfrog__" %in% colnames(d)) mean(d[, "n_leapfrog__"]) else NA_real_
  )
}


# diagnostic --------------------------------------------------------------

plot_batch_diagnostics <- function(diag_tbl) {
  p1 <- plot_ly(diag_tbl, x = ~ensemble, y = ~n_divergent, type = "bar", name = "divergences")
  p2 <- plot_ly(diag_tbl, x = ~ensemble, y = ~max_rhat, type = "scatter", mode = "markers", name = "max Rhat")
  p3 <- plot_ly(diag_tbl, x = ~ensemble, y = ~min_ess_bulk, type = "scatter", mode = "markers", name = "min ESS bulk")

  subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE) %>%
    layout(
      xaxis = list(title = "ensemble"),
      yaxis = list(title = "n_divergent"),
      yaxis2 = list(title = "max Rhat"),
      yaxis3 = list(title = "min ESS bulk")
    )
}

plot_sampler_cost <- function(diag_tbl) {
  plot_ly(diag_tbl, x = ~mean_n_leapfrog, y = ~mean_stepsize,
          type = "scatter", mode = "markers",
          text = ~paste("ensemble:", ensemble,
                        "<br>div:", n_divergent,
                        "<br>treedepth:", mean_treedepth)) %>%
    layout(
      xaxis = list(title = "mean n_leapfrog (cost driver)"),
      yaxis = list(title = "mean stepsize"),
      title = "Sampler cost vs stepsize"
    )
}


# batch runner ------------------------------------------------------------

run_batch <- function(rrs_hw_sub, mod, stan_data_base, out_dir,
                      mode = "sample",
                      chains = 2, iter_warmup = 300, iter_sampling = 300,
                      adapt_delta = 0.85, max_treedepth = 12,
                      threads_per_chain = 1) {

  rrs_hw_sub %>%
    mutate(
      stan = map2(
        data, ensemble,
        ~run_one_ensemble(
          df = .x,
          ensemble_id = .y,
          mod = mod,
          stan_data_base = stan_data_base,
          out_dir = out_dir,
          mode = mode,
          chains = chains,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          threads_per_chain = threads_per_chain
        )
      )
    )
}

run_batch_parallel <- function(rrs_hw_sub, mod, stan_data_base, out_dir,
                               mode = "sample",
                               chains = 2, iter_warmup = 300, iter_sampling = 300,
                               adapt_delta = 0.85, max_treedepth = 12,
                               threads_per_chain = 1,
                               workers = NULL,
                               seed = TRUE) {

  stopifnot(all(c("data", "ensemble") %in% names(rrs_hw_sub)))
  n <- nrow(rrs_hw_sub)
  if (n == 0) return(dplyr::mutate(rrs_hw_sub, stan = list()))

  if (!requireNamespace("future", quietly = TRUE))   stop("Install future")
  if (!requireNamespace("furrr", quietly = TRUE))    stop("Install furrr")
  if (!requireNamespace("progressr", quietly = TRUE)) stop("Install progressr")

  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  furrr_opts <- furrr::furrr_options(
    seed = seed,
    scheduling = 1,
    stdout = TRUE,
    #globals = TRUE # Ship function to worker via load_all
    packages = "SABER"
  )

  # ETA bookkeeping on the main R session (safe with progressr)
  times <- numeric(0)

  out <- progressr::with_progress({
    p <- progressr::progressor(steps = n)

    stan_list <- furrr::future_map2(
      rrs_hw_sub$data,
      rrs_hw_sub$ensemble,
      ~{
        t_i <- Sys.time()

        res <- run_one_ensemble(
          df = .x,
          ensemble_id = .y,
          mod = mod,
          stan_data_base = stan_data_base,
          out_dir = out_dir,
          mode = mode,
          chains = chains,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          threads_per_chain = threads_per_chain
        )

        dt <- as.numeric(difftime(Sys.time(), t_i, units = "secs"))

        # update ETA from completed tasks
        times <<- c(times, dt)
        done <- length(times)
        avg  <- mean(times)
        eta  <- (n - done) * avg

        p(sprintf("done %d/%d | last %.1fs | avg %.1fs/job | ETA ~%.0fs",
                  done, n, dt, avg, eta))

        res
      },
      .options = furrr_opts
    )

    rrs_hw_sub <- rrs_hw_sub %>% ungroup()

    dplyr::mutate(rrs_hw_sub, stan = stan_list)
  })

  out
}

extract_batch_tables <- function(results_tbl) {
  par_tbl <- results_tbl %>%
    transmute(ensemble, par_mean = map(stan, "par_mean")) %>%
    unnest(par_mean) %>%
    rename_with(~paste0(.x, "_estim"), .cols = -ensemble)

  diag_tbl <- results_tbl %>%
    transmute(
      ensemble,
      diag = map(stan, "diag")
    ) %>%
    unnest(diag)

  metric_tbl <- results_tbl %>%
    transmute(ensemble, metrics = map(stan, ~ .x$rrs$metrics)) %>%
    unnest(metrics)

  list(par = par_tbl, diag = diag_tbl, metrics = metric_tbl)
}

