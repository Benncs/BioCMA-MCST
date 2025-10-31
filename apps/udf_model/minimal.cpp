#include "Kokkos_Core.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_sampling/metropolis.hpp"
#include "mc/alias.hpp"
#include "mc/macros.hpp"
#include "mc/prng/prng_extension.hpp"
#include "models/udf_model.hpp"
#include "models/utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <mc/traits.hpp>
#include <udf_includes.hpp>

namespace
{
  using namespace Models;
  using FloatType = Models::UdfModel::FloatType;

  void __attribute__((constructor)) on_load()
  {
    std::printf("[UDF]: Minimal model loaded\r\n"); // NOLINT
  }

  constexpr FloatType glucose_to_biomass_yield = 0.5;

  constexpr FloatType l_dot_max = 2e-6 / 3600.; // m
  constexpr FloatType l_max_m = 2e-6;           // m
  constexpr FloatType l_min_m = l_max_m / 2.;   // m
  constexpr FloatType k = 1e-3;                 // m
  constexpr FloatType d_m = 0.6e-6;             // m
  constexpr FloatType lin_density =
      c_linear_density(static_cast<FloatType>(1000), d_m);

  constexpr FloatType phi_s_max =
      (l_dot_max * lin_density) / glucose_to_biomass_yield; // kg/s
  // Distributions
  const auto l_initial_dist = MC::Distributions::TruncatedNormal<float>(
      l_max_m * 0.75, l_max_m / 10., l_min_m, l_max_m);

  enum class particle_var : uint8_t
  {
    length = 0,
    l_max,
    phi_s,
    __COUNT__
  };

  std::size_t _set_nvar()
  {
    return static_cast<size_t>(particle_var::__COUNT__);
  };

  void _init_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                 std::size_t idx,
                 const Models::UdfModel::SelfParticle& arr,
                 const UdfModel::Config& config)
  {
    GET_PROPERTY(particle_var::length) = config(idx, 0);
    GET_PROPERTY(particle_var::l_max) = l_max_m;
  };

  MC::Status
  _update_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
              [[maybe_unused]] float d_t,
              std::size_t idx,
              const Models::UdfModel::SelfParticle& arr,
              const MC::LocalConcentration& c)
  {
    const auto s = static_cast<FloatType>(c(0));
    const FloatType g = s / (k + s);
    const FloatType phi_s = phi_s_max * g;
    const FloatType ldot = l_dot_max * g;
    const FloatType d_length = d_t * ldot;
    GET_PROPERTY(particle_var::length) += d_length / (1.0 + d_t * ldot);
    GET_PROPERTY(particle_var::phi_s) = -phi_s;

    return check_div(GET_PROPERTY(particle_var::length),
                     GET_PROPERTY(particle_var::l_max));
  }

  MC::ContribIndexBounds _get_bounds_udf()
  {
    int begin = INDEX_FROM_ENUM(particle_var::phi_s);
    return {.begin = begin, .end = begin + 1};
  }

  void
  _division_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                [[maybe_unused]] std::size_t idx,
                [[maybe_unused]] std::size_t idx2,
                [[maybe_unused]] const MC::DynParticlesModel<float>& arr,
                [[maybe_unused]] const MC::DynParticlesModel<float>& buffer_arr)
  {
    const FloatType new_current_length =
        GET_PROPERTY(particle_var::length) / 2.F;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length) =
        new_current_length;

    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = l_max_m;

    GET_PROPERTY(particle_var::length) = new_current_length;
    GET_PROPERTY(particle_var::l_max) = l_max_m;
  }

  double mass(std::size_t idx, const Models::UdfModel::SelfParticle& arr)
  {
    return GET_PROPERTY(particle_var::length) * lin_density;
  }

  std::vector<std::string_view> _names()
  {
    return {"length"};
  };

  std::vector<std::size_t> _get_number()
  {
    return {INDEX_FROM_ENUM(particle_var::length)};
  }

  UdfModel::Config _get_config_udf(const std::size_t n)
  {
    (void)n;
    Kokkos::ScopeGuard g{};

    UdfModel::Config res("config", n, 1);

    const auto target = KOKKOS_LAMBDA(const float x)
    {
      return (x < 2e-6 && x > 1e-6) ? 1. / 1e-6 : 0.;
    };

    Sampling::metropolis(
        target, Kokkos::subview(res, Kokkos::ALL, 0), 0.f, float(2e-6));

    return res;
  }
} // namespace
// clang-format off
EXPORT_MODULE(
    module, &_init_udf,
    &_update_udf,
    &_division_udf,
    &mass,
    &_names,
    &_get_number,
    &_set_nvar,
    &_get_bounds_udf,
    &_get_config_udf);
//clang-format on
