#include "saber_shim_stan.hpp"
#include <saber/saber.h>   // your C API

extern "C" double rt_model_c(const double* wl, int n_wl,
                            const double* a0, const double* a1,
                            const double* bb_w, int n_class,
                            const double* r_rs_b,
                            const double* /*params*/) {
  // 1) Build the spectral cache at these wavelengths
  build_cache(wl, (size_t)n_wl);

  // 2) As a minimal placeholder, compute the average bottom reflectance:
  const double* bottom = get_r_rs_b();
  double sum = 0.0;
  for (int i = 0; i < n_wl * n_class; ++i)
    sum += bottom[i];
  return sum / (n_wl * n_class);
}
