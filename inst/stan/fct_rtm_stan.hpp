#pragma once
#include <stan/math.hpp>
#include <Eigen/Dense>
#include "../include/rtm_core.hpp"


// Keep Stan-friendly signatures (return_type_t, pstream__) but call the core.
template <typename Tchl, typename Tag440, typename Tbbp550, typename Ts, typename Tgamma>
inline Eigen::Matrix<stan::return_type_t<Tchl,Tag440,Tbbp550,Ts,Tgamma>, Eigen::Dynamic, 2>
iop_from_oac_all(const Eigen::Ref<const Eigen::VectorXd>& wavelength,
                 const Eigen::Ref<const Eigen::VectorXd>& a_w,
                 const Eigen::Ref<const Eigen::VectorXd>& a0,
                 const Eigen::Ref<const Eigen::VectorXd>& a1,
                 const Eigen::Ref<const Eigen::VectorXd>& bb_w,
                 const Tchl& chl,
                 const Tag440& a_gnap_440,
                 const Tbbp550& bb_p_550,
                 const Ts& a_gnap_s,
                 const Tgamma& bb_p_gamma,
                 std::ostream* pstream__) {
  using T = stan::return_type_t<Tchl,Tag440,Tbbp550,Ts,Tgamma>;
  return saber::iop_from_oac_all_core<T>(wavelength, a_w, a0, a1, bb_w,
                                         T(chl), T(a_gnap_440), T(bb_p_550),
                                         T(a_gnap_s), T(bb_p_gamma));
}

template <typename DerivedA, typename DerivedBB, typename DerivedRB>
inline Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, 1>
forward_am03_ad(const Eigen::Ref<const Eigen::VectorXd>& wavelength,
                const Eigen::MatrixBase<DerivedA>& a,
                const Eigen::MatrixBase<DerivedBB>& bb,
                int water_type,
                double theta_sun_deg,
                double theta_view_deg,
                int shallow,
                const typename DerivedA::Scalar& h_w,
                const Eigen::MatrixBase<DerivedRB>& r_b,
                std::ostream* pstream__) {
  using T = typename DerivedA::Scalar;
  return saber::forward_am03_core<T>(wavelength,
                                     a.derived(), bb.derived(),
                                     water_type, theta_sun_deg, theta_view_deg,
                                     shallow, h_w, r_b.derived());
}
