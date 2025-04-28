#ifndef __MODELS_UPTAKES_HPP__
#define __MODELS_UPTAKES_HPP__

#include "mc/prng/prng_extension.hpp"
#include <cstddef>
#include <mc/traits.hpp>

#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

namespace Models
{
  enum Uptakeparticle_var : int
  {
    a_pts = 0,
    a_permease,
    n_permease,
    COUNT
  };

  // TODO Put ELSEWHERE
  template <FloatingPointType F> consteval F freq(F tau)
  {
    return F(1) / tau;
  }

  template <ModelType M> struct Uptake
  {

    static constexpr std::size_t n_var = static_cast<std::size_t>(Uptakeparticle_var::COUNT);
    static constexpr std::string_view name = "uptake";
    using uniform_weight = std::true_type;
    using Self = Uptake;
    using FloatType = M::FloatType;
    using SelfParticle = MC::ParticlesModel<M::n_var, Self::FloatType>;
    using Config = std::nullopt_t;

    MODEL_CONSTANT FloatType NPermease_init = 1;
    MODEL_CONSTANT FloatType NPermease_max = 200;
    MODEL_CONSTANT FloatType k_pts = 1e-3;
    MODEL_CONSTANT FloatType kppermease = 1e-2;

    MODEL_CONSTANT double tau_new_permease = 40.;
    MODEL_CONSTANT double tau_rm_perm = 200.;
    MODEL_CONSTANT double tau_pts = 20.;
    MODEL_CONSTANT double tau_Au = 40.;
    MODEL_CONSTANT double tau_Ad = 5.;

    static KOKKOS_INLINE_FUNCTION constexpr FloatType _k_pts()
    {
      return static_cast<FloatType>(k_pts);
    }

    static KOKKOS_INLINE_FUNCTION constexpr FloatType
    phi_pts(const FloatType phi_pts_max, const FloatType a_pts, const FloatType S)
    {
      return a_pts * (phi_pts_max * S / (S + static_cast<FloatType>(k_pts)));
    }

    static KOKKOS_INLINE_FUNCTION constexpr FloatType phi_permease(FloatType phi_permease_specific,
                                                                   FloatType n_permease,
                                                                   FloatType a_permease,
                                                                   FloatType S)
    {
      return a_permease * ((n_permease * phi_permease_specific) * (S) / ((S) + (kppermease)));
    }

    KOKKOS_INLINE_FUNCTION static void
    init(const MC::KPRNG::pool_type& random_pool, std::size_t idx, const SelfParticle& arr)
    {

      static constexpr FloatType half = FloatType(0.5);
      auto gen = random_pool.get_state();
      GET_PROPERTY(Uptakeparticle_var::a_pts) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen, FloatType(1e-3), FloatType(1e-4), FloatType(0.), FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::a_permease) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen, FloatType(0.8), FloatType(0.1), FloatType(0.), FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::n_permease) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(NPermease_init) * half,
              FloatType(NPermease_init) / FloatType(5.),
              FloatType(0.),
              FloatType(NPermease_init) * FloatType(10.));
      random_pool.free_state(gen);
    }

    KOKKOS_INLINE_FUNCTION static FloatType uptake(FloatType phi_pts_max,
                                                   FloatType phi_permease_specific,
                                                   std::size_t idx,
                                                   const SelfParticle& arr,
                                                   const MC::LocalConcentration& c)
    {
      const auto s = c(0);
      const FloatType phi_s_pts = phi_pts(phi_pts_max, GET_PROPERTY(Uptakeparticle_var::a_pts), s);
      return phi_s_pts + phi_permease(phi_permease_specific,
                                      GET_PROPERTY(Uptakeparticle_var::a_permease),
                                      GET_PROPERTY(Uptakeparticle_var::a_permease),
                                      s);
    }

    KOKKOS_INLINE_FUNCTION static FloatType uptake_step(FloatType phi_pts_max,
                                                        FloatType phi_permease_specific,
                                                        FloatType d_t,
                                                        std::size_t idx,
                                                        const SelfParticle& arr,
                                                        const MC::LocalConcentration& c)
    {

      constexpr FloatType apts_frequency = freq(tau_pts);
      constexpr FloatType aperm_frequency_au = freq(tau_Au);
      constexpr FloatType aperm_frequency_ad = freq(tau_Ad);
      constexpr FloatType n_perm_frequence_new = freq(tau_new_permease);
      constexpr FloatType n_perm_frequence_rm = freq(tau_rm_perm);

      const auto s = c(0);
      const FloatType phi_s_pts = phi_pts(phi_pts_max, GET_PROPERTY(Uptakeparticle_var::a_pts), s);
      const FloatType gamma_PTS_S = phi_s_pts / phi_pts_max;
      const FloatType phi_s = phi_s_pts + phi_permease(phi_permease_specific,
                                                       GET_PROPERTY(Uptakeparticle_var::n_permease),
                                                       GET_PROPERTY(Uptakeparticle_var::a_permease),
                                                       s);
      GET_PROPERTY(Uptakeparticle_var::a_pts) +=
          d_t * apts_frequency *
          (MONOD_RATIO(FloatType(1), s, k_pts) - GET_PROPERTY(Uptakeparticle_var::a_pts));

      GET_PROPERTY(Uptakeparticle_var::a_permease) +=
          d_t *
          ((aperm_frequency_au * gamma_PTS_S + aperm_frequency_ad * (FloatType(1) - gamma_PTS_S)) *
           (FloatType(1) - gamma_PTS_S - GET_PROPERTY(Uptakeparticle_var::a_permease)));

      GET_PROPERTY(Uptakeparticle_var::n_permease) +=
          d_t *
          (MONOD_RATIO(n_perm_frequence_new, k_pts, s) +
           MONOD_RATIO(n_perm_frequence_rm, s, k_pts)) *
          (MONOD_RATIO(NPermease_max, k_pts, s) - GET_PROPERTY(Uptakeparticle_var::n_permease));

      return phi_s;
    }

    // TODO FIND OUT HOW TO MAKE THIS MODEL
    //  KOKKOS_INLINE_FUNCTION static MC::Status update(const MC::KPRNG::pool_type& random_pool,
    //                                                  FloatType d_t,
    //                                                  std::size_t idx,
    //                                                  const SelfParticle& arr,
    //                                                  const MC::LocalConcentration& c)
    //  {
    //    constexpr FloatType apts_frequency = freq(tau_pts);
    //    constexpr FloatType aperm_frequency_au = freq(tau_Au);
    //    constexpr FloatType aperm_frequency_ad = freq(tau_Ad);
    //    constexpr FloatType n_perm_frequence_new = freq(tau_new_permease);
    //    constexpr FloatType n_perm_frequence_rm = freq(tau_rm_perm);
    //    const auto s = c(0);

    //   GET_PROPERTY(Uptakeparticle_var::a_pts) +=
    //       d_t * apts_frequency *
    //       (MONOD_RATIO(FloatType(1), s, k_pts) - GET_PROPERTY(Uptakeparticle_var::a_pts));
    //   GET_PROPERTY(Uptakeparticle_var::a_permease) +=
    //       d_t *
    //       ((aperm_frequency_au * gamma_PTS_S + aperm_frequency_ad * (FloatType(1) - gamma_PTS_S))
    //       *
    //        (FloatType(1) - gamma_PTS_S - GET_PROPERTY(Uptakeparticle_var::a_permease)));

    //   GET_PROPERTY(Uptakeparticle_var::n_permease) +=
    //       d_t *
    //       (MONOD_RATIO(n_perm_frequence_new, k_pts, s) +
    //        MONOD_RATIO(n_perm_frequence_rm, s, k_pts)) *
    //       (MONOD_RATIO(NPermease_init, k_pts, s) - GET_PROPERTY(Uptakeparticle_var::n_permease));
    // }

    KOKKOS_INLINE_FUNCTION static void division(const MC::KPRNG::pool_type& random_pool,
                                                std::size_t idx,
                                                std::size_t idx2,
                                                const SelfParticle& arr,
                                                const SelfParticle& buffer_arr)
    {

      static constexpr FloatType half = FloatType(0.5);
      auto generator = random_pool.get_state();
      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::a_pts) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              generator,
              GET_PROPERTY(Uptakeparticle_var::a_pts),
              GET_PROPERTY(Uptakeparticle_var::a_pts) * half,
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::a_permease) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              generator,
              GET_PROPERTY(Uptakeparticle_var::a_pts),
              GET_PROPERTY(Uptakeparticle_var::a_permease) * half,
              FloatType(0.),
              FloatType(1.));

      const auto new_n_permease = GET_PROPERTY(Uptakeparticle_var::n_permease) * half;
      GET_PROPERTY(Uptakeparticle_var::n_permease) = new_n_permease;
      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::n_permease) = new_n_permease;
      random_pool.free_state(generator);
    }

    [[maybe_unused]] KOKKOS_INLINE_FUNCTION static void
    contribution([[maybe_unused]] std::size_t idx,
                 [[maybe_unused]] std::size_t position,
                 [[maybe_unused]] double weight,
                 [[maybe_unused]] const SelfParticle& arr,
                 [[maybe_unused]] const MC::ContributionView& contributions)
    {
    }

    KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx, const SelfParticle& arr)
    {
      (void)idx;
      (void)arr;
      return 0.;
    }

    inline constexpr static std::array<std::string_view, n_var> names()
    {
      constexpr std::array<std::string_view, n_var> _names = {"a_pts", "a_permease", "n_permease"};
      static_assert(_names.size() == n_var);
      return _names;
    }
  };
  // CHECK_MODEL(Uptake<DefaultModel>)
} // namespace Models

#endif