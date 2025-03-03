#ifndef __BIOMODEL_UPTAKE__
#define __BIOMODEL_UPTAKE__

#include <Kokkos_Core.hpp>
#include <common/traits.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng_extension.hpp>
#include <type_traits>

#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

namespace Models::Uptake
{
#define DECLARE_UPTAKE_PROPERTIES(__type__)                                                        \
  __type__ a_pts;                                                                                  \
  __type__ a_permease;                                                                             \
  __type__ n_permease;

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

  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION constexpr F _k_pts()
  {
    return static_cast<F>(k_pts);
  }

  template <FloatingPointType F>
  KOKKOS_INLINE_FUNCTION F phi_pts(const F phi_pts_max, const F a_pts, const F S)
  {
    return a_pts * ((phi_pts_max) * (S) / ((S) + (_k_pts<F>())));
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

  template <FloatingPointType F, UptakeModel<F> Model, typename KokkosRNGPoolType>
  KOKKOS_INLINE_FUNCTION void
  distribute_division(Model& model, Model& child, KokkosRNGPoolType& generator)
  {
    static constexpr F half = F(0.5);
    // child.a_pts = Kokkos::min(1., Kokkos::max(generator.normal(model.a_pts, model.a_pts / 2.),
    // 0.));

    // child.a_permease =
    //     Kokkos::min(1., Kokkos::max(generator.normal(model.a_permease, model.a_permease / 2.),
    //     0.));

    child.a_pts = MC::Distributions::TruncatedNormal<F>::draw_from(
        generator, model.a_pts, model.a_pts * half, F(0.), F(1.));

    child.a_permease = MC::Distributions::TruncatedNormal<F>::draw_from(
        generator, model.a_pts, model.a_permease * half, F(0.), F(1.));

    model.n_permease = model.n_permease * half;
    child.n_permease = model.n_permease;
  }

  template <FloatingPointType F, UptakeModel<F> Model, typename KokkosRNGPoolType>
  KOKKOS_INLINE_FUNCTION void
  distribute_init(Model& newp, KokkosRNGPoolType& generator, F NPermease_init = 1)
  {
    static constexpr F half = F(0.5);
    // newp.a_permease = Kokkos::max(generator.normal(1e-3, 1e-4), 0.);
    // newp.a_pts = Kokkos::min(1., Kokkos::max(generator.normal(0.8, 0.1), 0.));
    // newp.n_permease = Kokkos::max(generator.normal(NPermease_init / 2., NPermease_init / 5.),
    // 0.);

    newp.a_permease =
        MC::Distributions::TruncatedNormal<F>::draw_from(generator, F(1e-3), F(1e-4), F(0.), F(1.));
    newp.a_pts =
        MC::Distributions::TruncatedNormal<F>::draw_from(generator, F(0.8), F(0.1), F(0.), F(1.));
    newp.n_permease = MC::Distributions::TruncatedNormal<F>::draw_from(generator,
                                                                       F(NPermease_init) * half,
                                                                       F(NPermease_init) / F(5.),
                                                                       F(0.),
                                                                       F(NPermease_init) * F(10.));
  }

} // namespace Models::Uptake

#endif