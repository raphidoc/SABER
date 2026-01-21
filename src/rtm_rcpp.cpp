// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "rtm_core.hpp"

// [[Rcpp::export]]
Eigen::MatrixXd iop_from_oac(const Eigen::VectorXd& wavelength,
                                   const Eigen::VectorXd& a_w,
                                   const Eigen::VectorXd& a0,
                                   const Eigen::VectorXd& a1,
                                   const Eigen::VectorXd& bb_w,
                                   double chl,
                                   double a_gnap_440,
                                   double bb_p_550,
                                   double a_gnap_s,
                                   double bb_p_gamma) {
  return saber::iop_from_oac_core<double>(wavelength, a_w, a0, a1, bb_w,
                                              chl, a_gnap_440, bb_p_550,
                                              a_gnap_s, bb_p_gamma);
}

//´ Albert & Mobley (2003) Forward Model
//´
//´ Computes Rrs below surface using analytical model from Albert & Mobley (2003)
//´
//´ @param wavelength vector of wavelengths [nm]
//´ @param iop list with elements `a` and `bb`, same length as wavelength
//´ @param water_type either 1 or 2 (default = 2)
//´ @param theta_sun sun zenith angle [degrees]
//´ @param theta_view sensor view angle [degrees]
//´ @param h_w optional: water depth [m] (enables shallow mode)
//´ @param r_b optional: bottom reflectance vector (same length as wavelength)
//´
//´ @return numeric vector of subsurface Rrs
//´
//´ @references Albert, A. and Mobley, C.D. (2003) ‘An analytical model for subsurface irradiance and remote sensing reflectance in deep and shallow case-2 waters’, Optics Express, 11(22), pp. 2873–2890. Available at: https://doi.org/10.1364/OE.11.002873.
// [[Rcpp::export]]
Eigen::VectorXd forward_am03(const Eigen::VectorXd& wavelength,
                               const Eigen::VectorXd& a,
                               const Eigen::VectorXd& bb,
                               int water_type,
                               double theta_sun_deg,
                               double theta_view_deg,
                               bool shallow,
                               double h_w,
                               Rcpp::Nullable<Eigen::VectorXd> r_b = R_NilValue) {
  const int n = wavelength.size();
  if (a.size() != n || bb.size() != n)
    Rcpp::stop("forward_am03_r: size mismatch (a/bb vs wavelength)");

  // If shallow=TRUE, need r_b; otherwise ignore it.
  Eigen::VectorXd rb_vec;
  if (shallow) {
    if (r_b.isNull())
      Rcpp::stop("forward_am03_r: shallow=TRUE requires r_b (length n_wl)");
    rb_vec = Rcpp::as<Eigen::VectorXd>(r_b);
    if (rb_vec.size() != n)
      Rcpp::stop("forward_am03_r: r_b size mismatch vs wavelength");
  } else {
    rb_vec.resize(0);
  }

  // Convert inputs to templated types expected by core
  Eigen::Matrix<double, Eigen::Dynamic, 1> a_t  = a;
  Eigen::Matrix<double, Eigen::Dynamic, 1> bb_t = bb;

  return saber::forward_am03_core<double>(
    wavelength,
    a_t, bb_t,
    water_type,
    theta_sun_deg,
    theta_view_deg,
    shallow ? 1 : 0,
    h_w,
    rb_vec
  );
}

// [[Rcpp::export]]
Eigen::VectorXd solve_rb_am03(const Eigen::VectorXd& wavelength,
                                const Eigen::VectorXd& a,
                                const Eigen::VectorXd& bb,
                                int water_type,
                                double theta_sun_deg,
                                double theta_view_deg,
                                double h_w,
                                const Eigen::VectorXd& rrs_obs) {
  const int n = wavelength.size();
  if (a.size() != n || bb.size() != n)
    Rcpp::stop("solve_rb_am03_r: size mismatch (a/bb vs wavelength)");
  if (rrs_obs.size() != n)
    Rcpp::stop("solve_rb_am03_r: rrs_obs size mismatch vs wavelength");

  // Core expects Eigen::Matrix<double,Dynamic,1> for a, bb
  Eigen::Matrix<double, Eigen::Dynamic, 1> a_t  = a;
  Eigen::Matrix<double, Eigen::Dynamic, 1> bb_t = bb;

  // rrs_obs can stay as Eigen::VectorXd; core reads rrs_obs(i) and casts to T
  return saber::solve_rb_am03_core<double>(
    wavelength,
    a_t,
    bb_t,
    water_type,
    theta_sun_deg,
    theta_view_deg,
    h_w,
    rrs_obs
  );
}
