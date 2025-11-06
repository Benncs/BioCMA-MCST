#ifndef __MODEL_MONOD_HPP__
#define __MODEL_MONOD_HPP__

#include "mc/macros.hpp"
#include "models/utils.hpp"
#include <Kokkos_MathematicalConstants.hpp>
#include <mc/alias.hpp>
#include <mc/traits.hpp>
#include <optional>
namespace Models
{

  struct Monod
  {
    enum class particle_var : int
    {
      l = 0,
      l_max,
      mu,
      _init_only_cell_lenghtening,
      phi_s_c,
      __COUNT__
    };

    static constexpr std::size_t n_var =
        INDEX_FROM_ENUM(particle_var::__COUNT__);
    static constexpr std::string_view name = "monod";
    using uniform_weight = std::true_type; // Using type alias
    using Self = Monod;
    using FloatType = float;
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
    using Config = std::nullopt_t;

    MODEL_CONSTANT FloatType y_s_x = 2;              // m
    MODEL_CONSTANT FloatType mu_max = 0.77 / 3600.;  // m
    MODEL_CONSTANT FloatType tau_meta = 1. / mu_max; // m
    MODEL_CONSTANT FloatType l_max_m = 2e-6;         // m
    MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.; // m
    MODEL_CONSTANT FloatType k_s = 1e-3;             // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;           // m
    MODEL_CONSTANT FloatType lin_density =
        c_linear_density(static_cast<FloatType>(1000), d_m);

    KOKKOS_INLINE_FUNCTION static void
    init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr)
    {
      auto generator = random_pool.get_state();
      // static constexpr auto mu_dist =
      // MC::Distributions::TruncatedNormal<float>(
      //     0.5/3600, 0.5/3600/5, 0., mu_max);

      static constexpr auto l_dist =
          MC::Distributions::TruncatedNormal<FloatType>(
              l_max_m * 0.75, l_max_m * 0.75 / 4, l_min_m, l_max_m);

      const auto l0 = l_dist.draw(generator);
      const auto mu_0 = mu_max; // mu_dist.draw(generator);
      random_pool.free_state(generator);

      GET_PROPERTY(particle_var::l) = l0;
      GET_PROPERTY(particle_var::l_max) = l_max_m;
      GET_PROPERTY(particle_var::mu) = mu_0;
      GET_PROPERTY(particle_var::_init_only_cell_lenghtening) =
          l0 / 2. / Kokkos::numbers::log2e_v<FloatType>;
    }

    KOKKOS_INLINE_FUNCTION static double
    mass([[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::l) * lin_density;
    }

    KOKKOS_INLINE_FUNCTION static MC::Status
    update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] FloatType d_t,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const MC::LocalConcentration& c)
    {
      const FloatType s = static_cast<FloatType>(Kokkos::max(0., c(0)));
      const FloatType mu_p = mu_max * s / (k_s + s);
      const FloatType mu_eff =
          Kokkos::min(GET_PROPERTY(Self::particle_var::mu), mu_p);

      GET_PROPERTY(Self::particle_var::l) +=
          d_t * (mu_eff *
                 GET_PROPERTY(Self::particle_var::_init_only_cell_lenghtening));
      GET_PROPERTY(Self::particle_var::mu) +=
          d_t * (1.0 / tau_meta) *
          (mu_p - GET_PROPERTY(Self::particle_var::mu));
      GET_PROPERTY(Self::particle_var::phi_s_c) =
          -mu_eff * y_s_x * mass(idx, arr);

      return check_div(GET_PROPERTY(Self::particle_var::l),
                       GET_PROPERTY(Self::particle_var::l_max));
    }

    KOKKOS_INLINE_FUNCTION static void
    division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
             [[maybe_unused]] std::size_t idx,
             [[maybe_unused]] std::size_t idx2,
             [[maybe_unused]] const SelfParticle& arr,
             [[maybe_unused]] const SelfParticle& buffer_arr)
    {
      const FloatType new_current_length = GET_PROPERTY(particle_var::l) / 2.F;
      GET_PROPERTY(particle_var::l) = new_current_length;

      // Newborn cell
      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l) = new_current_length;
      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = l_max_m;
      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::mu) =
          GET_PROPERTY(particle_var::mu);
      GET_PROPERTY_FROM(
          idx2, buffer_arr, particle_var::_init_only_cell_lenghtening) =
          GET_PROPERTY(particle_var::_init_only_cell_lenghtening);
    }

    static MC::ContribIndexBounds get_bounds()
    {
      int begin = INDEX_FROM_ENUM(Self::particle_var::phi_s_c);
      return {.begin = begin, .end = begin + 1};
    }

    static std::vector<std::string_view> names()
    {
      return {

          "length", "mu"};
    }

    static std::vector<std::size_t> get_number()
    {
      return {INDEX_FROM_ENUM(particle_var::l),
              INDEX_FROM_ENUM(particle_var::mu)};
    }
  };

  CHECK_MODEL(Monod)
} // namespace Models
#endif
