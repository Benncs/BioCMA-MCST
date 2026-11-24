#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include "mc/alias.hpp"
#include "mc/macros.hpp"
#include "models/udf_model.hpp"
#include "models/utils.hpp"
#include <common/env_var.hpp>
#include <cstdio>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <udf_includes.hpp>

/**
 * @brief Division time mode
 *
 * Two internal variables, age and length. Length grows from PTS glucose
 * uptake, division is on the age alone, so a starved cell still divides.
 * t_div is drawn once at birth from the chosen law, update only compares
 * against it.
 *
 *   BIOMC_UDF_DIVISION_LAW    0 deterministic (default), 1 exponential,
 *                             2 normal
 *   BIOMC_UDF_DIVISION_SIGMA  s, law 2 only, default tau_division / 10
 */
namespace
{
  using namespace Models;
  using FloatType = Models::UdfModel::FloatType;

  enum DivisionLaw : int
  {
    Deterministic = 0,
    Exponential,
    Normal
  };

  constexpr FloatType tau_division = 1200.; // s
  constexpr FloatType y_xs = 0.5;           // kgX/kgS
  constexpr FloatType k_pts = 1e-3;         // kg/m3
  constexpr FloatType l_min_m = 1e-6;       // m
  constexpr FloatType l_max_m = 2e-6;       // m
  constexpr FloatType d_m = 0.6e-6;         // m
  constexpr FloatType lin_density
      = c_linear_density(static_cast<FloatType>(1000), d_m);

  // A cell spans l_min to l_max over one interdivision time
  constexpr FloatType l_dot_max = (l_max_m - l_min_m) / tau_division;
  constexpr FloatType phi_s_pts_max = l_dot_max * lin_density / y_xs; // kgS/s

  // A draw below this would divide the cell on its first step
  constexpr FloatType t_div_floor = 1e-3; // s

  // Read once at load, host only, the udf backend is never cuda
  int division_law = DivisionLaw::Deterministic;
  FloatType division_sigma = tau_division / 10.;

  enum class particle_var : uint8_t
  {
    t_age = 0,
    t_div,
    length,
    phi_s,
    __COUNT__
  };

  void __attribute__((constructor))
  on_load()
  {
    division_law = Common::read_env_or("BIOMC_UDF_DIVISION_LAW", division_law);
    division_sigma
        = Common::read_env_or("BIOMC_UDF_DIVISION_SIGMA", division_sigma);

    std::printf("[UDF]: Division time model loaded, law=%d tau=%g sigma=%g\r\n",
                division_law,
                static_cast<double>(tau_division),
                static_cast<double>(division_sigma)); // NOLINT
  }

  /// @brief Inverse CDF of the chosen law, one draw per birth, never per step
  template <class Generator>
  KOKKOS_INLINE_FUNCTION FloatType
  draw_t_div(Generator& gen)
  {
    switch (division_law)
    {
    case DivisionLaw::Exponential:
      return Kokkos::max(
          MC::Distributions::Exponential<FloatType>{ 1.F / tau_division }.draw(
              gen),
          t_div_floor);
    case DivisionLaw::Normal:
      return Kokkos::max(
          MC::Distributions::Normal<FloatType>{ tau_division, division_sigma }
              .draw(gen),
          t_div_floor);
    default:
      return tau_division;
    }
  }

  std::size_t
  _set_nvar()
  {
    return static_cast<std::size_t>(particle_var::__COUNT__);
  };

  std::size_t
  _set_nc()
  {
    return 1;
  };

  void
  _init_udf(const MC::pool_type& random_pool,
            std::size_t idx,
            const Models::UdfModel::SelfParticle& arr,
            [[maybe_unused]] const UdfModel::Config& config)
  {
    auto gen = random_pool.get_state();
    GET_PROPERTY(particle_var::t_div) = draw_t_div(gen);
    random_pool.free_state(gen);

    GET_PROPERTY(particle_var::t_age) = 0.;
    GET_PROPERTY(particle_var::length) = l_min_m;
    GET_PROPERTY(particle_var::phi_s) = 0.;
  };

  MC::Status
  _update_udf([[maybe_unused]] const MC::pool_type& random_pool,
              float d_t,
              std::size_t idx,
              const Models::UdfModel::SelfParticle& arr,
              const Models::UdfModel::SelfContribs& arr_contribs,
              const std::size_t position_index,
              const MC::LocalConcentration& c)
  {
    const auto s = GET_CLAMPED_CONCENTRATION_CAST(FloatType, 0);
    const FloatType phi_s = phi_s_pts_max * s / (k_pts + s);
    GET_PROPERTY(particle_var::length) += d_t * phi_s * y_xs / lin_density;
    GET_PROPERTY(particle_var::t_age) += d_t;
    GET_PROPERTY(particle_var::phi_s) = phi_s;

    GET_CONTRIBS(0) = -phi_s;

    return check_div(GET_PROPERTY(particle_var::t_age),
                     GET_PROPERTY(particle_var::t_div));
  }

  void
  _division_udf(const MC::pool_type& random_pool,
                std::size_t idx,
                std::size_t idx2,
                const MC::DynParticlesModel<float>& arr,
                const MC::DynParticlesModel<float>& buffer_arr)
  {
    const FloatType new_current_length
        = GET_PROPERTY(particle_var::length) / 2.F;

    auto gen = random_pool.get_state();
    const FloatType t_div_1 = draw_t_div(gen);
    const FloatType t_div_2 = draw_t_div(gen);
    random_pool.free_state(gen);

    GET_PROPERTY(particle_var::length) = new_current_length;
    GET_PROPERTY(particle_var::t_age) = 0.;
    GET_PROPERTY(particle_var::t_div) = t_div_1;

    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length)
        = new_current_length;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::t_age) = 0.;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::t_div) = t_div_2;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::phi_s)
        = GET_PROPERTY(particle_var::phi_s);
  }

  double
  mass(std::size_t idx, const Models::UdfModel::SelfParticle& arr)
  {
    return GET_PROPERTY(particle_var::length) * lin_density;
  }

  std::vector<std::string_view>
  _names()
  {
    return { "t_age", "t_div", "length", "phi_s" };
  };

  std::vector<std::size_t>
  _get_number()
  {
    return { INDEX_FROM_ENUM(particle_var::t_age),
             INDEX_FROM_ENUM(particle_var::t_div),
             INDEX_FROM_ENUM(particle_var::length),
             INDEX_FROM_ENUM(particle_var::phi_s) };
  }

  std::vector<std::string_view>
  _species()
  {
    return { "S" };
  }

  UdfModel::Config
  _get_config_udf([[maybe_unused]] Kokkos::DefaultHostExecutionSpace& ep,
                  [[maybe_unused]] const std::size_t n)
  {
    return {};
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
    &_set_nc,
    &_get_config_udf,
    &_species);
//clang-format on
