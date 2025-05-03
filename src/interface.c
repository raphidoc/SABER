// src/interface.c
#include <R.h>
#include <Rinternals.h>
#include "saber.h"

SEXP c_load_pure_water(SEXP s_wl, SEXP s_a_w) {
  int n = LENGTH(s_wl);
  const double *wl = REAL(s_wl);
  const double *a_w = REAL(s_a_w);

  load_pure_water(wl, a_w, n);

  return R_NilValue;
}

SEXP c_load_a0_a1(SEXP s_wl, SEXP s_a0, SEXP s_a1) {
  int n = LENGTH(s_wl);
  const double *wl = REAL(s_wl);
  const double *a0 = REAL(s_a0);
  const double *a1 = REAL(s_a1);

  load_a0_a1(wl, a0, a1, n);

  return R_NilValue;
}

SEXP c_load_r_rs_b(SEXP s_wl,
                   SEXP s_class_names,
                   SEXP s_r_rs_b_matrix)
{
  /* 0.  Basic sanity checks ----------------------------------- */
  if (!isReal(s_wl) || !isString(s_class_names) || !isReal(s_r_rs_b_matrix))
    error("Invalid argument types.");

  R_xlen_t n_wl    = XLENGTH(s_wl);
  R_xlen_t n_class = XLENGTH(s_class_names);

  if (XLENGTH(s_r_rs_b_matrix) != n_wl * n_class)
    error("Matrix length (%lld) ≠ n_wl * n_class (%lld).",
          (long long)XLENGTH(s_r_rs_b_matrix),
          (long long)(n_wl * n_class));

  /* 1.  Borrow data from R ------------------------------------ */
  const double *wl     = REAL(s_wl);
  const double *matrix = REAL(s_r_rs_b_matrix);

  /* 2.  Build array of C‑string pointers (lives until function ends) */
  const char **class_ptrs =
  (const char **) R_alloc(n_class, sizeof(char*));  /* R‑managed */

  for (R_xlen_t i = 0; i < n_class; ++i)
    class_ptrs[i] = CHAR(STRING_ELT(s_class_names, i));  /* NO shadowing */

  /* 3.  Hand off to the C library (it deep‑copies the names) */
  if (load_r_rs_b(wl, class_ptrs, matrix,
                  (size_t)n_wl, (size_t)n_class) != 0)
    error("load_r_rs_b() failed (allocation error).");

  return R_NilValue;
}

SEXP c_iop_from_oac(SEXP s_wl, SEXP s_par) {
  int n = LENGTH(s_wl);
  const double *wl = REAL(s_wl);

  int n_param = LENGTH(s_par);
  const double *pval_ptr = REAL(s_par);
  SEXP names = getAttrib(s_par, R_NamesSymbol);

  const char **pnames = (const char **) R_alloc(n_param, sizeof(const char *));
  double *pvals = (double *) R_alloc(n_param, sizeof(double));

  for (int i = 0; i < n_param; i++) {
    pnames[i] = CHAR(STRING_ELT(names, i));
    pvals[i]  = pval_ptr[i];
  }

  double *a_out = (double *) R_alloc(n, sizeof(double));
  double *bb_out = (double *) R_alloc(n, sizeof(double));

  // Rprintf("Calling iop_from_oac with %d wavelengths and %d parameters\n", n, n_param);
  // for (int i = 0; i < n_param; i++) {
  //   Rprintf("  param[%d]: %s = %f\n", i, pnames[i], pval_ptr[i]);
  // }

  int status = iop_from_oac(wl, (size_t)n, pnames, pval_ptr, (size_t)n_param, a_out, bb_out);
  if (status != 0) {
    Rf_error("iop_from_oac failed (code %d).", status);
  }

  // Create R vectors
  SEXP a_vec = PROTECT(allocVector(REALSXP, n));
  SEXP bb_vec = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(a_vec), a_out, n * sizeof(double));
  memcpy(REAL(bb_vec), bb_out, n * sizeof(double));

  // Build named list
  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, a_vec);
  SET_VECTOR_ELT(result, 1, bb_vec);

  SEXP out_names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(out_names, 0, mkChar("a"));
  SET_STRING_ELT(out_names, 1, mkChar("bb"));
  setAttrib(result, R_NamesSymbol, out_names);

  UNPROTECT(4);
  return result;
}

SEXP c_compute_r_rs_b_lmm(SEXP s_fractions) {
  if (!isReal(s_fractions) || LENGTH(s_fractions) == 0) {
    Rf_error("`fractions` must be a non-empty numeric vector");
  }
  int n_cls = LENGTH(s_fractions);
  const double *fractions = REAL(s_fractions);

  // Extract class names from the names() attribute
  SEXP s_names = getAttrib(s_fractions, R_NamesSymbol);
  if (s_names == R_NilValue || LENGTH(s_names) != n_cls) {
    Rf_error("`fractions` must be a *named* numeric vector");
  }

  const char **class_names = (const char **) R_alloc(n_cls, sizeof(char *));
  for (int i = 0; i < n_cls; i++) {
    class_names[i] = CHAR(STRING_ELT(s_names, i));
  }

  int n_wl = get_n_wl(); // from your saber-lib

  SEXP result = PROTECT(allocVector(REALSXP, n_wl));
  double *rrs_b = REAL(result);

  compute_r_rs_b_lmm(class_names, fractions, (size_t)n_cls, rrs_b);

  UNPROTECT(1);
  return result;
}

SEXP c_forward_am03(SEXP s_wavelength, SEXP s_iop,
                    SEXP s_water_type,
                    SEXP s_theta_sun, SEXP s_theta_view,
                    SEXP s_h_w, SEXP s_r_b) {

  // ---- Wavelength and size ----
  int n = LENGTH(s_wavelength);
  const double *wl = REAL(s_wavelength);

  // ---- IOP: list with a and bb ----
  if (!isNewList(s_iop) || LENGTH(s_iop) != 2)
    Rf_error("`iop` must be a list with elements 'a' and 'bb'");

  SEXP names = getAttrib(s_iop, R_NamesSymbol);
  if (names == R_NilValue || LENGTH(names) != 2)
    Rf_error("`iop` list must be named with 'a' and 'bb'");

  SEXP s_a = R_NilValue, s_bb = R_NilValue;
  for (int i = 0; i < 2; i++) {
    const char *name = CHAR(STRING_ELT(names, i));
    if (strcmp(name, "a") == 0) s_a = VECTOR_ELT(s_iop, i);
    else if (strcmp(name, "bb") == 0) s_bb = VECTOR_ELT(s_iop, i);
  }
  if (s_a == R_NilValue || s_bb == R_NilValue)
    Rf_error("`iop` must contain both 'a' and 'bb'");

  const double *a = REAL(s_a);
  const double *bb = REAL(s_bb);

  // ---- Angles and water type ----
  int water_type = INTEGER(s_water_type)[0];
  double theta_sun = REAL(s_theta_sun)[0];
  double theta_view = REAL(s_theta_view)[0];

  // ---- Depth (optional) ----
  double h = (s_h_w == R_NilValue) ? 0.0 : REAL(s_h_w)[0];

  // ---- Bottom reflectance (optional) ----
  int shallow = (s_r_b != R_NilValue);
  const double *r_b = shallow ? REAL(s_r_b) : NULL;

  // ---- Output vector ----
  SEXP result = PROTECT(allocVector(REALSXP, n));
  double *rrs = REAL(result);

  forward_am03(
    wl, a, bb, (size_t)n,
    water_type, theta_sun, theta_view,
    shallow, h, r_b, rrs
  );

  UNPROTECT(1);
  return result;
}

