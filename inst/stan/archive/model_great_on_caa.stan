functions {
  // Returns [n_wl, 2] where col 1 = a(λ), col 2 = bb(λ)
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

  // Pre-interpolated LUTs on the same wavelength grid
  vector[n_wl] a_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
  vector[n_wl] bb_w;

  // Bottom library on the same wavelength grid
  int<lower=1> n_class;
  matrix[n_wl, n_class] r_b_lib;
  array[3] int<lower=1, upper=n_class> bottom_class_ids;

  // Geometry / flags
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

  real<lower=0, upper=15> h_w;

  //simplex[3] bottom_mix; // composition (sums to 1)
  //real<lower=0> a_bottom; // overall bottom reflectance amplitude

  vector<lower=0>[3] bottom_w;   // free nonnegative weights

  //real<lower=0> z_gate;   // depth where g=0.5
  //real<lower=0.5> z_width;  // softness (m)

  real<lower=0> sigma_model; // model discrepancy (same units as Rrs)
}

transformed parameters {
  vector[n_wl] a;
  vector[n_wl] bb;
  vector[n_wl] r_b;
  vector[n_wl] rrs_model;
  vector[n_wl] sigma_tot;

  // IOPs from OAC (autodiff-able)
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

  // Bottom reflectance mixture with amplitude scale
  //vector[n_wl] r_b_shape;

  //r_b_shape =
      //bottom_mix[1] * to_vector(r_b_lib[, bottom_class_ids[1]]) +
      //bottom_mix[2] * to_vector(r_b_lib[, bottom_class_ids[2]]) +
      //bottom_mix[3] * to_vector(r_b_lib[, bottom_class_ids[3]]);

  // optional but recommended: normalize shape so A_bottom is interpretable
  //r_b_shape = r_b_shape / mean(r_b_shape);

  //r_b = a_bottom * r_b_shape;

  //
  r_b =
      bottom_w[1] * to_vector(r_rs_b_lib[, bottom_class_ids[1]]) +
      bottom_w[2] * to_vector(r_rs_b_lib[, bottom_class_ids[2]]) +
      bottom_w[3] * to_vector(r_rs_b_lib[, bottom_class_ids[3]]);

  // Test depth gating bottom reflectance contribution
  //real g = inv_logit( (z_gate - h_w) / z_width );

  //vector[n_wl] rrs_deep = forward_am03_ad(
  //  wavelength, a, bb,
  //  water_type, theta_sun_deg, theta_view_deg,
  //  0, h_w, rep_vector(0.0, n_wl) // r_b ignored in deep mode anyway
  //);

  //vector[n_wl] rrs_shallow = forward_am03_ad(
  //  wavelength, a, bb,
  //  water_type, theta_sun_deg, theta_view_deg,
  //  1, h_w, r_b
  //);

  //rrs_model = (1 - g) * rrs_deep + g * rrs_shallow;


  // Forward RT model
  rrs_model = forward_am03_ad(
    wavelength, a, bb,
    water_type, theta_sun_deg, theta_view_deg,
    shallow, h_w, r_b
  );

  // total uncertainty per band
  sigma_tot = sqrt(square(sigma_rrs) + square(sigma_model));
}

model {
  // Priors
  chl       ~ lognormal(log(0.1), 0.7);
  a_g_440   ~ lognormal(log(0.5), 0.7);
  a_nap_440 ~ lognormal(log(0.5), 0.7);
  bb_p_550  ~ lognormal(log(0.005), 0.7);

  a_g_s      ~ lognormal(log(0.017), 0.25);
  a_nap_s    ~ lognormal(log(0.0116), 0.25);
  bb_p_gamma ~ lognormal(log(0.46), 0.25);

  h_w ~ normal(5, 3);

  // Bottom reflectance magnitude uncertainty in mean reflectance unit.
  //a_bottom ~ lognormal(log(0.18), 0.25);
  //a_bottom ~ normal(0.018, 0.02);

  // Optional: weakly-informative prior favoring not-too-extreme mixtures
  // bottom_mix ~ dirichlet(rep_vector(1.0, 3));

  bottom_w ~ normal(0.18, 0.4);

  // Depth gate prior
  //z_gate  ~ normal(7.0, 2.0);   // bottom matters at mu, sigma
  //z_width ~ lognormal(log(1.0), 0.3);   // transition scale; keep >0 via <lower=0>

  // Pick a scale related to typical Rrs magnitudes in your dataset.
  sigma_model ~ normal(0, 0.00005); // half-normal due to <lower=0>

  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_tot);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
}
