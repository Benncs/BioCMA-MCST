#ifndef __MODELS_UPTAKES_DYN_HPP__
#define __MODELS_UPTAKES_DYN_HPP__

#include "Kokkos_Assert.hpp"
#include <common/traits.hpp>
#include <cstddef>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/uptake.hpp>

enum Uptakeparticle_var : int
{
  a_pts = 0,
  a_permease_1,
  a_permease_2,
  COUNT
};

// TODO Put ELSEWHERE
template <FloatingPointType F> consteval F freq(F tau)
{
  return F(1) / tau;
}

namespace Models
{

  // TODO Put ELSEWHERE

  template <FloatingPointType F>
  KOKKOS_INLINE_FUNCTION constexpr F f_saturation(const F x, const F k) noexcept
  {
    return x / (x + k);
  }

  template <typename T>
  concept UptakeModel = requires(T model) {
    // { T::NPermease_init } -> std::convertible_to<typename T::FloatType>;

    // { T::NPermease_max } -> std::convertible_to<typename T::FloatType>;

    { T::k } -> std::convertible_to<typename T::FloatType>;
    { T::k_perm } -> std::convertible_to<typename T::FloatType>;

    { T::tau_new_permease } -> std::convertible_to<typename T::FloatType>;

    { T::tau_rm_perm } -> std::convertible_to<typename T::FloatType>;

    { T::tau_pts } -> std::convertible_to<typename T::FloatType>;

    { T::tau_Au } -> std::convertible_to<typename T::FloatType>;

    { T::tau_Ad } -> std::convertible_to<typename T::FloatType>;

    { T::beta } -> std::convertible_to<typename T::FloatType>;

    //    { T::delta } -> std::convertible_to<typename T::FloatType>;
  };

  template <UptakeModel U, ModelType M = U> struct Uptake
  {
    static constexpr std::size_t n_var =
        static_cast<std::size_t>(Uptakeparticle_var::COUNT);
    static constexpr std::string_view name = "uptake";
    using uniform_weight = std::true_type;
    using Self = Uptake;
    using FloatType = M::FloatType;
    using SelfParticle = MC::ParticlesModel<M::n_var, Self::FloatType>;
    using Config = std::nullopt_t;

    static KOKKOS_INLINE_FUNCTION constexpr FloatType
    f_G(const FloatType x) noexcept
    {
      return f_saturation(x, U::k);
    }

    static KOKKOS_INLINE_FUNCTION constexpr FloatType phi_pts(
        const FloatType phi_pts_max, const FloatType a_pts, const FloatType S)
    {
      return a_pts * phi_pts_max * f_G(S);
    }

    static KOKKOS_INLINE_FUNCTION constexpr FloatType
    phi_permease(FloatType phi_pts_max,
                 FloatType a_permease_2,
                 FloatType a_permease,
                 FloatType S)
    {
      constexpr FloatType k_perm = U::k / 10;
      constexpr FloatType sigma = k_perm / 10.;
      const FloatType innerexp = (S - k_perm) / sigma;
      return a_permease * a_permease_2 * phi_pts_max * U::beta *
             Kokkos::exp(-innerexp * innerexp);
    }

    KOKKOS_INLINE_FUNCTION static void
    init(const MC::KPRNG::pool_type& random_pool,
         std::size_t idx,
         const SelfParticle& arr)
    {

      // static constexpr FloatType half = FloatType(0.5);
      auto gen = random_pool.get_state();
      // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      GET_PROPERTY(Uptakeparticle_var::a_pts) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(1e-3),
              FloatType(1e-4),
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::a_permease_1) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(0.8),
              FloatType(0.1),
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::a_permease_2) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(0.8),
              FloatType(0.1),
              FloatType(0.),
              FloatType(1.));

      // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      random_pool.free_state(gen);
    }

    KOKKOS_INLINE_FUNCTION static FloatType
    uptake(FloatType phi_pts_max,
           std::size_t idx,
           const SelfParticle& arr,
           const MC::LocalConcentration& c,
           FloatType* r_phi_pts = nullptr,
           FloatType* r_phi_perm = nullptr)
    {
      const auto s = c(0);

      const FloatType phi_s_pts =
          phi_pts(phi_pts_max, GET_PROPERTY(Uptakeparticle_var::a_pts), s);
      if (r_phi_pts != nullptr)
      {
        *r_phi_pts = phi_s_pts;
      }
      const auto phi_s_perm =
          phi_permease(phi_pts_max,
                       GET_PROPERTY(Uptakeparticle_var::a_permease_2),
                       GET_PROPERTY(Uptakeparticle_var::a_permease_1),
                       s);
      if (r_phi_perm != nullptr)
      {
        *r_phi_perm = phi_s_perm;
      }

      return phi_s_pts + phi_s_perm;
    }

    KOKKOS_INLINE_FUNCTION static FloatType
    uptake_step(FloatType phi_pts_max,
                FloatType d_t,
                std::size_t idx,
                const SelfParticle& arr,
                const MC::LocalConcentration& c,
                FloatType* r_phi_pts = nullptr,
                FloatType* r_phi_perm = nullptr)
    {
      constexpr FloatType apts_frequency = 1. / U::tau_pts;
      constexpr FloatType aperm_frequency_au = 1. / U::tau_Au;
      constexpr FloatType aperm_frequency_ad = 1. / U::tau_Ad;
      constexpr FloatType n_perm_frequence_new = 1. / U::tau_new_permease;
      constexpr FloatType n_perm_frequence_rm = 1. / U::tau_rm_perm;
      constexpr FloatType one = FloatType(1);

      auto& a_pts = GET_PROPERTY(Uptakeparticle_var::a_pts);
      auto& a_p1 = GET_PROPERTY(Uptakeparticle_var::a_permease_1);
      auto& a_p2 = GET_PROPERTY(Uptakeparticle_var::a_permease_2);

      const auto s = c(0);
      KOKKOS_ASSERT(s >= 0.F);
      const auto G = f_G(s);
      const auto a_pts_g = a_pts * G;
      const auto phi_s_pts = phi_pts(phi_pts_max, a_pts, s);
      const auto phi_perm = phi_permease(phi_pts_max, a_p2, a_p1, s);
      const FloatType phi_s = phi_s_pts + phi_perm;
      if (r_phi_perm != nullptr)
      {
        *r_phi_perm = phi_perm;
      }
      if (r_phi_pts != nullptr)
      {
        *r_phi_pts = phi_s_pts;
      }

      a_pts += d_t * apts_frequency * (G - a_pts);

      a_p1 += d_t * ((aperm_frequency_au * a_pts_g +
                      aperm_frequency_ad * (one - a_pts_g)) *
                     (one - G - a_p1));

      a_p2 +=
          d_t * ((n_perm_frequence_rm * G + n_perm_frequence_new * (one - G)) *
                 (one - a_pts_g - a_p2));

      return phi_s;
    }

    KOKKOS_INLINE_FUNCTION static void
    division(const MC::KPRNG::pool_type& random_pool,
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

      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::a_permease_1) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              generator,
              GET_PROPERTY(Uptakeparticle_var::a_pts),
              GET_PROPERTY(Uptakeparticle_var::a_permease_1) * half,
              FloatType(0.),
              FloatType(1.));

      const auto new_n_permease =
          GET_PROPERTY(Uptakeparticle_var::a_permease_2) * half;
      GET_PROPERTY(Uptakeparticle_var::a_permease_2) = new_n_permease;
      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::a_permease_2) =
          new_n_permease;
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

    KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx,
                                              const SelfParticle& arr)
    {
      (void)idx;
      (void)arr;
      return 0.;
    }

    inline constexpr static std::array<std::string_view, n_var> names()
    {
      constexpr std::array<std::string_view, n_var> _names = {
          "a_pts", "a_permease_1", "a_permease_2"};
      static_assert(_names.size() == n_var);
      return _names;
    }
  };

  // CHECK_MODEL(Uptake<DefaultModel>)
} // namespace Models

#endif
