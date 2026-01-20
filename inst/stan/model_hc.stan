functions {
  // Returns [n_wl, 2] where col 1 = a(λ), col 2 = bb(λ)
  matrix iop_from_oac_all(
      vector wavelength,
      vector a_w,
      vector a0,
      vector a1,
      vector bb_w,
      real chl,
      real a_gnap_440,
      real bb_p_550,
      real a_gnap_s,
      real bb_p_gamma
  );

  vector forward_am03_ad(
      vector wavelength,
      vector a,
      vector bb,
      int water_type,
      real theta_sun_deg,
      real theta_view_deg,
      int shallow,
      real h_w,
      vector r_b
  );
}

data {
  int<lower=1> n_wl;
  vector[n_wl] wavelength;
  vector[n_wl] rrs_obs;
  vector<lower=1e-8>[n_wl] sigma_rrs;

  // Pre-interpolated LUTs on the same wavelength grid
  vector[n_wl] a_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
  vector[n_wl] bb_w;

  // Bottom library on the same wavelength grid
  int<lower=1> n_class;
  matrix[n_wl, n_class] r_b_mu_lib;
  matrix[n_wl, n_class] r_b_sd_lib;
  array[3] int<lower=1, upper=n_class> bottom_class_ids;

  // Geometry / flags
  int<lower=1,upper=2> water_type;
  real theta_sun;
  real theta_view;
  int<lower=0,upper=1> shallow;

  // KNOWN depth (passed from data)
  real<lower=0> h_w;
}

parameters {
  real<lower=0> chl;
  real<lower=0> a_gnap_440;
  real<lower=0> bb_p_550;

  real<lower=0> a_gnap_s;
  real<lower=0> bb_p_gamma;

  // real<lower=0> a_nap_s;
  // real<lower=0> a_g_s;
  // real<lower=0> bb_p_gamma;

  simplex[3] r_b_mix;
  real<lower=0> r_b_a;
  // vector[K] delta_raw;

  real<lower=0> sigma_model; // model uncertainty (same units as Rrs)
}

transformed parameters {
  vector[n_wl] a;
  vector[n_wl] bb;

  vector[n_wl] r_b_mu_mix;
  vector[n_wl] r_b_sd_mix;
  vector[n_wl] delta;

  vector[n_wl] r_b;
  vector[n_wl] rrs_model;
  vector[n_wl] sigma_tot;

  // IOPs from OAC (autodiff-able)
  {
    matrix[n_wl, 2] iop = iop_from_oac_all(
      wavelength, a_w, a0, a1, bb_w,
      chl, a_gnap_440, bb_p_550,
      a_gnap_s, bb_p_gamma
    );
    a  = iop[, 1];
    bb = iop[, 2];
  }

  // mixture mean
  r_b_mu_mix =
      r_b_mix[1] * to_vector(r_b_mu_lib[, bottom_class_ids[1]]) +
      r_b_mix[2] * to_vector(r_b_mu_lib[, bottom_class_ids[2]]) +
      r_b_mix[3] * to_vector(r_b_mu_lib[, bottom_class_ids[3]]);

  // mixture SD (variance mixture approx: sum w^2 * sd^2)
  //r_b_sd_mix = sqrt(
      //square(r_b_mix[1]) * square(to_vector(r_b_sd_lib[, bottom_class_ids[1]])) +
      //square(r_b_mix[2]) * square(to_vector(r_b_sd_lib[, bottom_class_ids[2]])) +
     //square(r_b_mix[3]) * square(to_vector(r_b_sd_lib[, bottom_class_ids[3]]))
    //);

  // deviation in SD units, low-rank
  // delta = Phi * (delta_scale * delta_raw);

  // bottom spectrum
  // r_b = r_b_a * (r_b_mu_mix + r_b_sd_mix .* delta);
  r_b = r_b_a * r_b_mu_mix;

  // Forward RT model
  rrs_model = forward_am03_ad(
    wavelength, a, bb,
    water_type, theta_sun, theta_view,
    shallow, h_w, r_b
  );

  // total uncertainty per band
  sigma_tot = sqrt(square(sigma_rrs) + square(sigma_model));
}

model {
  // Priors
  chl ~ lognormal(log(0.9), 0.6);
  a_gnap_440 ~ lognormal(log(0.32), 0.30);
  bb_p_550 ~ lognormal(log(0.0012), 0.45);

  a_gnap_s ~ normal(0.0168, 0.0012);
  bb_p_gamma ~ normal(0.45, 0.08);

  // bottom priors
  r_b_a     ~ lognormal(log(1.0), 0.5);
  // delta_raw ~ normal(0, 1);

  // model error
  sigma_model ~ normal(0, 0.00005); // half-normal due to <lower=0>

  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_tot);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
}

