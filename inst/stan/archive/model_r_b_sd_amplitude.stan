functions {
  matrix iop_from_oac_all(
      vector wavelength,
      vector a_w,
      vector a0,
      vector a1,
      vector bb_w,
      real chl,
      real a_g_440,
      real a_nap_440,
      real a_g_s,
      real a_nap_s,
      real bb_p_550,
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

  vector[n_wl] a_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
  vector[n_wl] bb_w;

  // Bottom library mean + sd on same wavelength grid
  int<lower=1> n_class;
  matrix[n_wl, n_class] r_b_mu_lib;
  matrix<lower=0>[n_wl, n_class] r_b_sd_lib;
  array[3] int<lower=1, upper=n_class> bottom_class_ids;

  int<lower=1,upper=2> water_type;
  real theta_sun_deg;
  real theta_view_deg;
  int<lower=0,upper=1> shallow;
}

parameters {
  real<lower=0> chl;
  real<lower=0, upper=10> a_g_440;
  real<lower=0> a_nap_440;
  real<lower=0> bb_p_550;

  real<lower=0> a_nap_s;
  real<lower=0> a_g_s;
  real<lower=0> bb_p_gamma;

  real<lower=0, upper=30> h_w;

  simplex[3] bottom_mix;
  real<lower=0> a_bottom;

  // latent bottom reflectance spectrum used by forward model
  vector<lower=0, upper=1>[n_wl] r_b_latent;

  real<lower=0> sigma_model;
}

transformed parameters {
  vector[n_wl] a;
  vector[n_wl] bb;

  // mixture mean & sd (in reflectance units, before amplitude scaling)
  vector[n_wl] r_b_mu_shape;
  vector[n_wl] r_b_sd_shape;

  // scaled versions (after amplitude)
  vector[n_wl] r_b_mu_mix;
  vector[n_wl] r_b_sd_mix;

  vector[n_wl] rrs_model;
  vector[n_wl] sigma_tot;

  // IOPs
  {
    matrix[n_wl, 2] iop = iop_from_oac_all(
      wavelength, a_w, a0, a1, bb_w,
      chl, a_g_440, a_nap_440,
      a_g_s, a_nap_s,
      bb_p_550, bb_p_gamma
    );
    a  = iop[, 1];
    bb = iop[, 2];
  }

  // Bottom mixture mean
  r_b_mu_shape =
      bottom_mix[1] * to_vector(r_b_mu_lib[, bottom_class_ids[1]]) +
      bottom_mix[2] * to_vector(r_b_mu_lib[, bottom_class_ids[2]]) +
      bottom_mix[3] * to_vector(r_b_mu_lib[, bottom_class_ids[3]]);

  // Bottom mixture sd
  r_b_sd_shape =
      bottom_mix[1] * to_vector(r_b_sd_lib[, bottom_class_ids[1]]) +
      bottom_mix[2] * to_vector(r_b_sd_lib[, bottom_class_ids[2]]) +
      bottom_mix[3] * to_vector(r_b_sd_lib[, bottom_class_ids[3]]);

  // Optional normalization of SHAPE.
  // IMPORTANT: if you normalize the mean, you should scale sd accordingly.
  //{
    //real mu_scale = mean(r_b_mu_shape);
    //r_b_mu_shape = r_b_mu_shape / mu_scale;
    //r_b_sd_shape = r_b_sd_shape / mu_scale;
  //}

  // Apply amplitude scaling
  r_b_mu_mix = a_bottom * r_b_mu_shape;
  r_b_sd_mix = a_bottom * r_b_sd_shape;

  // Forward RT model uses latent bottom spectrum (not the mean)
  rrs_model = forward_am03_ad(
    wavelength, a, bb,
    water_type, theta_sun_deg, theta_view_deg,
    shallow, h_w, r_b_latent
  );

  sigma_tot = sqrt(square(sigma_rrs) + square(sigma_model));
}

model {
  // Priors (yours)
  chl       ~ lognormal(log(0.1), 0.7);
  a_g_440   ~ lognormal(log(0.5), 0.7);
  a_nap_440 ~ lognormal(log(0.5), 0.7);
  bb_p_550  ~ lognormal(log(0.005), 0.7);

  a_g_s      ~ lognormal(log(0.017), 0.25);
  a_nap_s    ~ lognormal(log(0.0116), 0.25);
  bb_p_gamma ~ lognormal(log(0.46), 0.25);

  // bottom magnitude prior (yours)
  //a_bottom ~ lognormal(log(0.0155), 0.15);

  // composition prior (I recommend turning it on; mild)
  //bottom_mix ~ dirichlet(rep_vector(2.0, 3));

  // model discrepancy (you may later make this data-scaled)
  sigma_model ~ normal(0, 0.0005);

  // ---- SOFT CONSTRAINT ON BOTTOM SPECTRUM ----
  // r_b_latent is constrained around the mixture mean with SD from measured uncertainty.
  // Add a small floor to avoid sd=0 causing numerical issues.
  r_b_latent ~ normal(r_b_mu_mix, r_b_sd_mix + 1e-6);

  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_tot);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
  vector[n_wl] r_b_hat = r_b_latent;
}
