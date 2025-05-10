#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP c_saber_version();

extern SEXP c_load_pure_water(SEXP, SEXP);
extern SEXP c_load_a0_a1(SEXP, SEXP, SEXP);
extern SEXP c_load_r_rs_b(SEXP, SEXP, SEXP);
extern SEXP c_iop_from_oac(SEXP, SEXP);
extern SEXP c_compute_r_rs_b_lmm(SEXP);
extern SEXP c_forward_am03(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_retrieve_r_rs_b_am03(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP c_build_cache(SEXP);
extern SEXP c_saber_reset_tables();
extern SEXP c_get_a_w();

static const R_CallMethodDef CallEntries[] = {
  {"c_saber_version",  (DL_FUNC) &c_saber_version,  0},
  {"c_load_pure_water", (DL_FUNC) &c_load_pure_water, 2},
  {"c_load_a0_a1", (DL_FUNC) &c_load_a0_a1, 3},
  {"c_load_r_rs_b", (DL_FUNC) &c_load_r_rs_b, 3},
  {"c_iop_from_oac", (DL_FUNC) &c_iop_from_oac, 2},
  {"c_compute_r_rs_b_lmm", (DL_FUNC) &c_compute_r_rs_b_lmm, 1},
  {"c_forward_am03", (DL_FUNC) &c_forward_am03, 7},
  {"c_retrieve_r_rs_b_am03", (DL_FUNC) &c_retrieve_r_rs_b_am03, 7},
  {"c_build_cache", (DL_FUNC) &c_build_cache, 1},
  {"c_saber_reset_tables", (DL_FUNC) &c_saber_reset_tables, 0},
  {"c_get_a_w", (DL_FUNC) &c_get_a_w, 0},
  {NULL, NULL, 0}
};

void R_init_SABER(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
