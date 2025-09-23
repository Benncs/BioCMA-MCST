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
  ap_1 = 0,
  ap_2,
  ap_3,
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

    { T::beta } -> std::convertible_to<typename T::FloatType>;
    { T::tau_ap_1 } -> std::convertible_to<typename T::FloatType>;
    { T::tau_ap_2 } -> std::convertible_to<typename T::FloatType>;
    { T::tau_ap_2 } -> std::convertible_to<typename T::FloatType>;

    // { T::tau_new_permease } -> std::convertible_to<typename T::FloatType>;

    // { T::tau_rm_perm } -> std::convertible_to<typename T::FloatType>;

    // { T::tau_pts } -> std::convertible_to<typename T::FloatType>;

    // { T::tau_Au } -> std::convertible_to<typename T::FloatType>;

    // { T::tau_Ad } -> std::convertible_to<typename T::FloatType>;

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

    static KOKKOS_INLINE_FUNCTION constexpr FloatType
    phi_pts(const FloatType phi_max, const FloatType a_pts, const FloatType S)
    {
      return a_pts * phi_max * f_G(S);
    }

    static KOKKOS_INLINE_FUNCTION constexpr FloatType
    phi_permease(FloatType phi_max, FloatType ap_2, FloatType ap_3, FloatType S)
    {
      return ap_2 * ap_3 * phi_max * U::beta * f_saturation(S, U::k_perm);
    }

    // static KOKKOS_INLINE_FUNCTION consexpr FloatType phi(FloatType phi_max,
    // FloatType s)
    //{
    // constexpr FloatType k_perm = U::k / 10;
    // constexpr FloatType sigma = k_perm / 10.;
    // const FloatType innerexp = (s - k_perm) / sigma;
    // return a_permease * a_permease_2 * phi_max * U::beta *
    // Kokkos::exp(-innerexp * innerexp);
    //
    // return phi_max * (a_pts * f_G(s) + )
    //}
    //
    KOKKOS_INLINE_FUNCTION
    static void init(const MC::KPRNG::pool_type& random_pool,
                     std::size_t idx,
                     const SelfParticle& arr)
    {

      // static constexpr FloatType half = FloatType(0.5);
      auto gen = random_pool.get_state();
      // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      GET_PROPERTY(Uptakeparticle_var::ap_1) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(1e-3),
              FloatType(1e-4),
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::ap_2) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              gen,
              FloatType(0.8),
              FloatType(0.1),
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY(Uptakeparticle_var::ap_3) =
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
    uptake(FloatType phi_max,
           std::size_t idx,
           const SelfParticle& arr,
           const MC::LocalConcentration& c,
           FloatType* r_phi_pts = nullptr,
           FloatType* r_phi_perm = nullptr)
    {
      const auto s = c(0);

      const FloatType phi_s_pts =
          phi_pts(phi_max, GET_PROPERTY(Uptakeparticle_var::ap_1), s);
      if (r_phi_pts != nullptr)
      {
        *r_phi_pts = phi_s_pts;
      }
      const auto phi_s_perm =
          phi_permease(phi_max,
                       GET_PROPERTY(Uptakeparticle_var::ap_1),
                       GET_PROPERTY(Uptakeparticle_var::ap_2),
                       s);
      if (r_phi_perm != nullptr)
      {
        *r_phi_perm = phi_s_perm;
      }

      return phi_s_pts + phi_s_perm;
    }

    KOKKOS_INLINE_FUNCTION static FloatType
    uptake_step(FloatType phi_max,
                FloatType d_t,
                std::size_t idx,
                const SelfParticle& arr,
                const MC::LocalConcentration& c,
                FloatType* r_phi_pts = nullptr,
                FloatType* r_phi_perm = nullptr)
    {
      constexpr FloatType one = FloatType(1);
      constexpr FloatType f_ap_1 = 1. / U::tau_ap_1;
      constexpr FloatType f_ap_2 = 1. / U::tau_ap_2;
      constexpr FloatType f_ap_3 = 1. / U::tau_ap_3;

      auto& ap_1 = GET_PROPERTY(Uptakeparticle_var::ap_1);
      auto& ap_2 = GET_PROPERTY(Uptakeparticle_var::ap_2);
      auto& ap_3 = GET_PROPERTY(Uptakeparticle_var::ap_3);

      auto phi_s = uptake(phi_max, idx, arr, c, r_phi_pts, r_phi_perm);

      const auto s = c(0);
      KOKKOS_ASSERT(s >= 0.F);
      const auto G = f_G(s);

      // constexpr FloatType k_perm = U::k_perm;
      // constexpr FloatType sigma = k_perm;
      // const FloatType innerexp = (s - 10*k_perm) / sigma;

      const FloatType ap1_star = G;
      //      const FloatType ap2_star = Kokkos::exp(-innerexp*innerexp);

      constexpr auto yon = U::k_perm / 10.;
      constexpr auto yoff = U::k_perm * 5.;
      constexpr auto alpha = (1. - 0.) / (yoff - yon);
      constexpr auto beta = -yon * alpha;

      const FloatType ap2_star = (s >= yon) & (s < yoff);
      const FloatType ap3_star = (one - ap_1 * G);

      ap_1 += d_t * f_ap_1 * (ap1_star - ap_1);

      ap_2 += d_t * f_ap_2 * (ap2_star - ap_2);

      ap_3 += d_t * f_ap_3 * (ap3_star - ap_3);

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
      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::ap_1) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              generator,
              GET_PROPERTY(Uptakeparticle_var::ap_1),
              GET_PROPERTY(Uptakeparticle_var::ap_1) * half,
              FloatType(0.),
              FloatType(1.));

      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::ap_2) =
          MC::Distributions::TruncatedNormal<FloatType>::draw_from(
              generator,
              GET_PROPERTY(Uptakeparticle_var::ap_1),
              GET_PROPERTY(Uptakeparticle_var::ap_2) * half,
              FloatType(0.),
              FloatType(1.));

      const auto new_n_permease = GET_PROPERTY(Uptakeparticle_var::ap_3) * half;
      GET_PROPERTY(Uptakeparticle_var::ap_3) = new_n_permease;
      GET_PROPERTY_FROM(idx2, buffer_arr, Uptakeparticle_var::ap_3) =
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
