#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

extern SEXP c_snell_law(SEXP, SEXP);

// Helper: convert degrees to radians
static inline double deg2rad(double deg) {
  return deg * M_PI / 180.0;
}

// [[register]]
SEXP c_forward_am03(SEXP wavelength_sexp, SEXP iop_list, SEXP water_type_sexp,
                    SEXP theta_sun_sexp, SEXP theta_view_sexp,
                    SEXP h_w_sexp, SEXP r_b_sexp) {

  int nprotect = 0;

  if (!Rf_isReal(wavelength_sexp)) Rf_error("wavelength must be numeric");
  if (!Rf_isNewList(iop_list)) Rf_error("iop must be a list with elements 'a' and 'bb'");
  if (!Rf_isReal(theta_sun_sexp) || !Rf_isReal(theta_view_sexp))
    Rf_error("theta_sun and theta_view must be numeric scalars");

  int n = LENGTH(wavelength_sexp);
  double* wl = REAL(wavelength_sexp);

  SEXP a_vec = VECTOR_ELT(iop_list, 0);
  SEXP bb_vec = VECTOR_ELT(iop_list, 1);
  if (LENGTH(a_vec) != n || LENGTH(bb_vec) != n)
    Rf_error("iop$a and iop$bb must be same length as wavelength");

  double* a = REAL(a_vec);
  double* bb = REAL(bb_vec);

  int water_type = INTEGER(water_type_sexp)[0]; // 1 or 2
  double theta_sun = deg2rad(REAL(theta_sun_sexp)[0]);
  double theta_view = deg2rad(REAL(theta_view_sexp)[0]);

  // Optional shallow parameters
  int shallow = (!Rf_isNull(h_w_sexp) && !Rf_isNull(r_b_sexp));
  double* r_b = NULL;
  double h_w = 0.0;

  if (shallow) {
    if (!Rf_isReal(h_w_sexp) || !Rf_isReal(r_b_sexp))
      Rf_error("h_w and r_b must be numeric if provided");
    if (LENGTH(r_b_sexp) != n)
      Rf_error("r_b must match length of wavelength");

    r_b = REAL(r_b_sexp);
    h_w = REAL(h_w_sexp)[0];
  }

  // Output vector
  SEXP rrs_out = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  double* rrs = REAL(rrs_out);

  SEXP geometry = PROTECT(c_snell_law(theta_view_sexp, theta_sun_sexp)); nprotect++;
  double view_w = REAL(VECTOR_ELT(geometry, 0))[0];
  double sun_w  = REAL(VECTOR_ELT(geometry, 1))[0];

  for (int i = 0; i < n; i++) {
    double ext = a[i] + bb[i];
    double omega_b = bb[i] / ext;

    // Fresnel & geometric effects (simplified snell model)
    double f_rs;
    if (water_type == 1) {
      f_rs = 0.095;
    } else {
      f_rs = 0.0512 *
        (1 + (4.6659 * omega_b) +
        (-7.8387 * omega_b * omega_b) +
        (5.4571 * omega_b * omega_b * omega_b)) *
        (1 + (0.1098 / cos(sun_w))) *
        (1 + (0.4021 / cos(view_w)));
    }

    double rrs_deep = f_rs * omega_b;

    if (shallow) {
      double k0 = (water_type == 1) ? 1.0395 : 1.0546;

      double Kd = k0 * (ext / cos(sun_w));
      double kuW = (ext / cos(view_w)) *
        pow(1 + omega_b, 3.5421) *
        (1 - (0.2786 / cos(sun_w)));

      double kuB = (ext / cos(view_w)) *
        pow(1 + omega_b, 2.2658) *
        (1 - (0.0577 / cos(sun_w)));

      double Ars1 = 1.1576;
      double Ars2 = 1.0389;

      rrs[i] = rrs_deep * (1 - (Ars1 * exp(-h_w * (Kd + kuW)))) +
        Ars2 * r_b[i] * exp(-h_w * (Kd + kuB));
    } else {
      rrs[i] = rrs_deep;
    }
  }

  UNPROTECT(nprotect);
  return rrs_out;
}
