functions {
  // Returns [n_wl, 2] where col 1 = a(λ), col 2 = bb(λ)
  matrix iop_from_oac_spm(
    vector wavelength,
    vector a_w,
    vector a0,
    vector a1,
    vector bb_w,
    real chl,
    real a_g_440,
    real spm,
    real a_nap_star,
    real bb_p_star,
    real a_g_s,
    real a_nap_s,
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

  // LUTs on same wavelength grid
  vector[n_wl] a_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
  vector[n_wl] bb_w;

  // ---- Bottom spline prior: class-conditional mixture on spline coefficients ----
  int<lower=1> K;
  int<lower=2> m;                 // number of spline basis functions
  matrix[n_wl, m] rb_B;           // spline basis evaluated on wavelength grid (include intercept!)

  array[K] vector[m] rb_c_mu_k;   // class means in coefficient space (logit space)
  array[K] matrix[m, m] rb_c_L_k; // lower Cholesky of coef covariance per class
  simplex[K] rb_pi;               // mixture weights

  // Geometry / flags
  int<lower=1,upper=2> water_type;
  real theta_sun;
  real theta_view;
  int<lower=0,upper=1> shallow;

  // prior hyperparameters
  real a_nap_star_mu;
  real a_nap_star_sd;
  real bb_p_star_mu;
  real bb_p_star_sd;

  real a_g_s_mu;
  real a_g_s_sd;
  real a_nap_s_mu;
  real a_nap_s_sd;
  real bb_p_gamma_mu;
  real bb_p_gamma_sd;

  // Depth prior
  real h_w_mu;
  real h_w_sd;
}

parameters {
  // Water / IOP parameters
  real<lower=0> chl;
  real<lower=0> a_g_440;
  real<lower=0> spm;

  real<lower=0> a_nap_star;
  real<lower=0> bb_p_star;

  real<lower=0> a_g_s;
  real<lower=0> a_nap_s;
  real<lower=0> bb_p_gamma;

  real<lower=0, upper=30> h_w;

  // Bottom: spline coefficients (logit space)
  vector[m] rb_c;

  // Model discrepancy (same units as Rrs)
  real<lower=0> sigma_model;
}

transformed parameters {
  vector[n_wl] a;
  vector[n_wl] bb;

  vector[n_wl] rb_eta;     // logit-space bottom spectrum
  vector[n_wl] r_b;        // bottom reflectance (0,1)

  vector[n_wl] rrs_model;
  vector[n_wl] sigma_tot;

  // IOPs (autodiff)
  {
    matrix[n_wl, 6] iop = iop_from_oac_spm(
      wavelength, a_w, a0, a1, bb_w,
      chl, a_g_440, spm,
      a_nap_star, bb_p_star,
      a_g_s, a_nap_s, bb_p_gamma
    );
    a  = iop[, 1];
    bb = iop[, 2];
  }

  // Bottom spectrum from spline coefficients
  rb_eta = rb_B * rb_c;
  r_b = inv_logit(rb_eta);

  // Forward model
  rrs_model = forward_am03_ad(
    wavelength, a, bb,
    water_type, theta_sun, theta_view,
    shallow, h_w, r_b
  );

  sigma_tot = sqrt(square(sigma_rrs) + square(sigma_model));
}

model {
  // ----- Priors: water / IOP -----
  chl      ~ lognormal(log(0.9), 0.6);
  a_g_440  ~ lognormal(log(0.25), 0.6);
  spm      ~ lognormal(log(2), 1.6);

  a_nap_star  ~ normal(a_nap_star_mu, a_nap_star_sd);
  bb_p_star   ~ normal(bb_p_star_mu, bb_p_star_sd);

  a_g_s       ~ normal(a_g_s_mu, a_g_s_sd);
  a_nap_s     ~ normal(a_nap_s_mu, a_nap_s_sd);
  bb_p_gamma  ~ normal(bb_p_gamma_mu, bb_p_gamma_sd);

  h_w ~ normal(h_w_mu, h_w_sd);

  // ---- Bottom prior: mixture on coefficients (shape+magnitude linked by class) ----
  {
    vector[K] log_comp;
    for (k in 1:K) {
      // rb_c ~ MVN(mu_k, L_k L_k')
      vector[m] z = mdivide_left_tri_low(rb_c_L_k[k], rb_c - rb_c_mu_k[k]);
      real lp = -0.5 * (m * log(2*pi())
                        + 2 * sum(log(diagonal(rb_c_L_k[k])))
                        + dot_self(z));
      log_comp[k] = log(rb_pi[k]) + lp;
    }
    target += log_sum_exp(log_comp);
  }

  // Discrepancy
  sigma_model ~ normal(0, 0.00005);

  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_tot);
}

generated quantities {
  vector[n_wl] r_b_hat = r_b;
  vector[n_wl] rrs_hat = rrs_model;
}
