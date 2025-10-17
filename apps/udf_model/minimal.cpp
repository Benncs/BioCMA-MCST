#include "mc/alias.hpp"
#include "mc/macros.hpp"
#include "mc/prng/prng_extension.hpp"
#include "models/utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <mc/traits.hpp>
#include <udf_includes.hpp>
using namespace Models;
using FloatType = Models::UdfModel::FloatType;

void __attribute__((constructor)) on_load()
{
  std::printf("[UDF]: Minimal model loaded\r\n"); // NOLINT
}

MODEL_CONSTANT FloatType glucose_to_biomass_yield = 0.5;

MODEL_CONSTANT FloatType l_dot_max = 2e-6 / 3600.; // m
MODEL_CONSTANT FloatType l_max_m = 2e-6;           // m
MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.;   // m
MODEL_CONSTANT FloatType k = 1e-3;                 // m
MODEL_CONSTANT FloatType d_m = 0.6e-6;             // m
MODEL_CONSTANT FloatType lin_density =
    c_linear_density(static_cast<FloatType>(1000), d_m);

MODEL_CONSTANT FloatType phi_s_max =
    (l_dot_max * lin_density) / glucose_to_biomass_yield; // kg/s
// Distributions
MODEL_CONSTANT auto l_initial_dist = MC::Distributions::TruncatedNormal<float>(
    l_min_m * 0.75, l_min_m * 0.75 / 2., l_min_m, l_max_m);

enum class particle_var : int
{
  age = 0,
  length,
  contrib_phi_s,
  __COUNT__
};

std::size_t _set_nvar()
{
  return static_cast<size_t>(particle_var::__COUNT__);
};

void _init_udf(const MC::KPRNG::pool_type& random_pool,
               std::size_t idx,
               const Models::UdfModel::SelfParticle& arr)
{

  auto gen = random_pool.get_state();
  const FloatType linit = l_initial_dist.draw(gen);
  random_pool.free_state(gen);
  GET_PROPERTY(particle_var::age) = 0.;
  GET_PROPERTY(particle_var::length) = linit;
};

MC::Status
_update_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
            [[maybe_unused]] float d_t,
            [[maybe_unused]] std::size_t idx,
            [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
            [[maybe_unused]] const MC::LocalConcentration& c)
{
  const auto s = static_cast<FloatType>(c(0));
  const FloatType g = s / (k + s);
  FloatType phi_s = phi_s_max * g;
  const auto ldot = l_dot_max * g;

  GET_PROPERTY(particle_var::age) += d_t;
  GET_PROPERTY(particle_var::length) += d_t * ldot;
  GET_PROPERTY(particle_var::contrib_phi_s) = phi_s;

  return (GET_PROPERTY(particle_var::length) >= l_max_m) ? MC::Status::Division
                                                         : MC::Status::Idle;
}

void _contribution_udf(
    [[maybe_unused]] std::size_t idx,
    [[maybe_unused]] std::size_t position,
    [[maybe_unused]] double weight,
    [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
    [[maybe_unused]] const MC::ContributionView& contributions)
{
  auto access = contributions.access();
  access(position, 0) +=
      -weight * GET_PROPERTY(particle_var::contrib_phi_s); // NOLINT
}

void _division_udf(
    [[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
    [[maybe_unused]] std::size_t idx,
    [[maybe_unused]] std::size_t idx2,
    [[maybe_unused]] const MC::DynParticlesModel<float>& arr,
    [[maybe_unused]] const MC::DynParticlesModel<float>& buffer_arr)
{
  const FloatType new_current_length = GET_PROPERTY(particle_var::length) / 2.F;
  GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length) =
      new_current_length;
  GET_PROPERTY(particle_var::length) = new_current_length;
}

double mass(std::size_t idx, const Models::UdfModel::SelfParticle& arr)
{
  return GET_PROPERTY(particle_var::length) * lin_density;
}

std::vector<std::string_view> _names()
{
  return {"age", "length"};
};

std::vector<std::size_t> _get_number()
{
  return {INDEX_FROM_ENUM(particle_var::age),
          INDEX_FROM_ENUM(particle_var::length)};
}

// clang-format off
EXPORT_MODULE(
    module, &_init_udf,
    &_update_udf,
    &_contribution_udf,
    &_division_udf,
    &mass,
    &_names,
    &_get_number,
    &_set_nvar);
//clang-format on
