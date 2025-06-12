#ifndef __MODELS_UPTAKES_DYN_HPP__
#define __MODELS_UPTAKES_DYN_HPP__

#include "common/common.hpp"
#include <common/traits.hpp>
#include <cstddef>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/uptake.hpp>
#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

namespace Models
{

  // TODO Put ELSEWHERE

  template <typename T>
  concept DynUptakeModel = ModelType<T> && requires(T model) {

    T::NPermease_init;
    T::NPermease_max ;
    T::k_pts ;
    T::kppermease ;
    T::tau_new_permease ;
    T::tau_rm_perm ;
    T::tau_pts ;
    T::tau_Au;
    T::tau_Ad;
  };


template <DynUptakeModel M> struct UptakeDyn
{

  static constexpr std::size_t n_var = static_cast<std::size_t>(Uptakeparticle_var::COUNT);
  static constexpr std::string_view name = "uptake";
  using uniform_weight = std::true_type;
  using Self = UptakeDyn;
  using FloatType = M::FloatType;
  using SelfParticle = MC::ParticlesModel<M::n_var, Self::FloatType>;
  using Config = std::nullopt_t;

  static KOKKOS_INLINE_FUNCTION constexpr FloatType
  phi_pts(const FloatType phi_pts_max, const FloatType a_pts, const FloatType S)
  {
    return a_pts * (phi_pts_max * S / (S + static_cast<FloatType>(M::k_pts)));
  }

  static KOKKOS_INLINE_FUNCTION constexpr FloatType phi_permease(FloatType phi_permease_specific,
                                                                 FloatType n_permease,
                                                                 FloatType a_permease,
                                                                 FloatType S)
  {
    return a_permease *
           ((n_permease * phi_permease_specific) * (S) / ((S) + (M::kppermease)));
  }

  KOKKOS_INLINE_FUNCTION static void init(const MC::KPRNG::pool_type& random_pool,
                                          std::size_t idx,
                                          const SelfParticle& arr)
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
            FloatType(M::NPermease_init) * half,
            FloatType(M::NPermease_init) / FloatType(5.),
            FloatType(0.),
            FloatType(M::NPermease_init) * FloatType(10.));
    random_pool.free_state(gen);
  }

  KOKKOS_INLINE_FUNCTION static FloatType uptake(FloatType phi_pts_max,
                                                 FloatType phi_permease_specific,
                                                 std::size_t idx,
                                                 const SelfParticle& arr,
                                                 const MC::LocalConcentration& c,
                                                 FloatType* r_phi_pts = nullptr,
                                                 FloatType* r_phi_perm = nullptr)
  {
    const auto s = c(0);

    const FloatType phi_s_pts = phi_pts(phi_pts_max, GET_PROPERTY(Uptakeparticle_var::a_pts), s);
    if (r_phi_pts != nullptr)
    {
      *r_phi_pts = phi_s_pts;
    }
    const auto phi_s_perm = phi_permease(phi_permease_specific,
                                         GET_PROPERTY(Uptakeparticle_var::a_permease),
                                         GET_PROPERTY(Uptakeparticle_var::a_permease),
                                         s);
    if (r_phi_perm != nullptr)
    {
      *r_phi_perm = phi_s_perm;
    }

    return phi_s_pts + phi_s_perm;
  }

  KOKKOS_INLINE_FUNCTION static FloatType uptake_step(FloatType phi_pts_max,
                                                      FloatType phi_permease_specific,
                                                      FloatType d_t,
                                                      std::size_t idx,
                                                      const SelfParticle& arr,
                                                      const MC::LocalConcentration& c,
                                                      FloatType* r_phi_pts = nullptr,
                                                      FloatType* r_phi_perm = nullptr)
  {
    constexpr FloatType apts_frequency = 1. / M::tau_pts;
    constexpr FloatType aperm_frequency_au = 1. / M::tau_Au;
    constexpr FloatType aperm_frequency_ad = 1. / M::tau_Ad;
    constexpr FloatType n_perm_frequence_new = 1. / M::tau_new_permease;
    constexpr FloatType n_perm_frequence_rm = 1. / M::tau_rm_perm;

    const auto s = c(0);
    const FloatType phi_s_pts = phi_pts(phi_pts_max, GET_PROPERTY(Uptakeparticle_var::a_pts), s);
    const FloatType gamma_PTS_S = phi_s_pts / phi_pts_max;

    const auto G = MONOD_RATIO(FloatType(1), s, M::k_pts);

    if (r_phi_pts != nullptr)
    {
      *r_phi_pts = phi_s_pts;
    }
    const auto phi_perm = phi_permease(phi_permease_specific,
                                       GET_PROPERTY(Uptakeparticle_var::n_permease),
                                       GET_PROPERTY(Uptakeparticle_var::a_permease),
                                       s);
    if (r_phi_perm != nullptr)
    {
      *r_phi_perm = phi_perm;
    }

    const FloatType phi_s = phi_s_pts + phi_perm;

    GET_PROPERTY(Uptakeparticle_var::a_pts) +=
        d_t * apts_frequency *
        (G - GET_PROPERTY(Uptakeparticle_var::a_pts));

    GET_PROPERTY(Uptakeparticle_var::a_permease) +=
        d_t *
        ((aperm_frequency_au * gamma_PTS_S + aperm_frequency_ad * (FloatType(1) - gamma_PTS_S)) *
         (FloatType(1) - G - GET_PROPERTY(Uptakeparticle_var::a_permease)));

    // GET_PROPERTY(Uptakeparticle_var::n_permease) +=
    //     d_t *
    //     (MONOD_RATIO(n_perm_frequence_new, M::k_pts, s) +
    //      MONOD_RATIO(n_perm_frequence_rm, s, M::k_pts)) *
    //     (MONOD_RATIO(M::NPermease_max, M::k_pts, s) -
    //      GET_PROPERTY(Uptakeparticle_var::n_permease));

    GET_PROPERTY(Uptakeparticle_var::n_permease) +=
        d_t *
        (MONOD_RATIO(n_perm_frequence_new, M::k_pts, s) +
         MONOD_RATIO(n_perm_frequence_rm, s, M::k_pts)) *
         (M::NPermease_max*(FloatType(1) - gamma_PTS_S) - GET_PROPERTY(Uptakeparticle_var::n_permease));

    return phi_s;
  }

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