#ifndef __BIOMODEL_UPTAKE__
#define __BIOMODEL_UPTAKE__

#include "Kokkos_Macros.hpp"
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

#define DECLARE_UPTAKE_PROPERTIES(__type__)                                                        \
  __type__ a_pts;                                                                                  \
  __type__ a_permease;                                                                             \
  __type__ n_permease;

namespace Models::Uptake
{

  template <typename T>
  concept FloatingPointType = std::is_floating_point_v<std::remove_reference_t<T>>;

  //   template <typename T, typename F>
  //   concept UptakeModel = FloatingPointType<F> && requires(T& obj) {
  //     { obj.a_pts } -> std::same_as<std::add_lvalue_reference_t<std::remove_reference_t<F>>>;
  //     { obj.n_permease } ->
  //     std::same_as<std::add_lvalue_reference_t<std::remove_reference_t<F>>>; { obj.a_permease }
  //     -> std::same_as<std::add_lvalue_reference_t<std::remove_reference_t<F>>>;
  //   };

  template <typename T, typename F>
  concept UptakeModel = FloatingPointType<F> && requires(T& obj) {
    { obj.a_pts } -> std::same_as<std::add_lvalue_reference_t<F>>;
    { obj.n_permease } -> std::same_as<std::add_lvalue_reference_t<F>>;
    { obj.a_permease } -> std::same_as<std::add_lvalue_reference_t<F>>;
  };

  struct Tau
  {
    static constexpr double new_permease = 40.;
    static constexpr double tau_rm_perm = 200.;
    static constexpr double pts = 20.;
    static constexpr double Au = 40.;
    static constexpr double Ad = 5.;
  };

  constexpr double k_pts = 1e-3;
  constexpr double kppermease = 1e-2;

  template <FloatingPointType F>
  KOKKOS_INLINE_FUNCTION F phi_pts(const F phi_pts_max, const F a_pts, const F S)
  {
    return a_pts * ((phi_pts_max) * (S) / ((S) + (k_pts)));
  }

  template <FloatingPointType F>
  KOKKOS_INLINE_FUNCTION F phi_permease(F phi_permease_specific, F n_permease, F a_permease, F S)
  {
    return a_permease * ((n_permease * phi_permease_specific) * (S) / ((S) + (kppermease)));
  }

  template <FloatingPointType F, UptakeModel<F> Model>
  KOKKOS_INLINE_FUNCTION void step(F d_t, F gamma_PTS_S, F s, Model& model, F NPermease_max)
  {
    model.a_pts += d_t * 1.0 / Tau::pts * (MONOD_RATIO(1., s, k_pts) - model.a_pts);
    model.a_permease +=
        d_t * (((1.0 / Tau::Au) * gamma_PTS_S + (1.0 / Tau::Ad) * (1.0 - gamma_PTS_S)) *
               (1.0 - gamma_PTS_S - model.a_permease));

    model.n_permease += d_t *
                        (MONOD_RATIO(1. / Tau::new_permease, k_pts, s) +
                         MONOD_RATIO(1. / Tau::tau_rm_perm, s, k_pts)) *
                        (MONOD_RATIO(NPermease_max, k_pts, s) - model.n_permease);
  }

  template <FloatingPointType F, UptakeModel<F> Model>
  KOKKOS_FUNCTION F
  uptake(F d_t, Model& model, F s, F phi_pts_max, F phi_permease_specific, F NPermease_max)
  {
    const F phi_s_pts = phi_pts(phi_pts_max, model.a_pts, s);
    const F gamma_PTS_S = phi_s_pts / phi_pts_max;
    const F phi_s =
        phi_s_pts + phi_permease(phi_permease_specific, model.n_permease, model.a_permease, s);

    step(d_t, gamma_PTS_S, s, model, NPermease_max);

    return phi_s;
  }


} // namespace Models::Uptake

#endif