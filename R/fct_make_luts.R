interp1 <- function(x, y, xout) {
  approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
}

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
