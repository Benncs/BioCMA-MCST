#include "mc/alias.hpp"
#include "mc/macros.hpp"
#include "mc/prng/prng_extension.hpp"
#include "models/utils.hpp"
#include <cstdlib>
#include <mc/traits.hpp>
#include <udf_includes.hpp>
using namespace Models;

void __attribute__((constructor)) on_load()
{
  printf("UDF loaded\n"); // NOLINT
}
using FloatType = Models::UdfModel::FloatType;
MODEL_CONSTANT FloatType l_max_m = 5e-6;   // m
MODEL_CONSTANT FloatType l_c_m = 3e-6;     // m
MODEL_CONSTANT FloatType d_m = 0.6e-6;     // m
MODEL_CONSTANT FloatType l_min_m = 0.9e-6; // m
MODEL_CONSTANT FloatType lin_density = c_linear_density(static_cast<FloatType>(1000), d_m);

enum class particle_var : int
{
  age = 0,
  length,
  l_cp,
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
  constexpr auto local_lc = 3e-6;
  constexpr auto length_dist_c = MC::Distributions::TruncatedNormal<Models::UdfModel::FloatType>(
      local_lc , l_c_m / 20., l_min_m, l_max_m);

  constexpr auto length_dist =
      MC::Distributions::Normal<Models::UdfModel::FloatType>(local_lc, 0.1e-6);

  auto gen = random_pool.get_state();
  GET_PROPERTY(particle_var::length) = 2.5e-6; // Kokkos::max(FloatType(0), length_dist.draw(gen));
  random_pool.free_state(gen);
  GET_PROPERTY(particle_var::l_cp) = length_dist_c.draw(gen);
  GET_PROPERTY(particle_var::contrib_phi_s) = 0;
};

MC::Status _update_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                       [[maybe_unused]] float d_t,
                       [[maybe_unused]] std::size_t idx,
                       [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
                       [[maybe_unused]] const MC::LocalConcentration& c)
{

  (void)random_pool;
  const auto phi_s = c(0);

  auto nu = 6e-19 * phi_s;

  auto dl = (phi_s > 1e-6) ? 1e-9 : 0;

  GET_PROPERTY(particle_var::contrib_phi_s) = -2. * nu*1e-5;

  GET_PROPERTY(particle_var::length) += static_cast<float>(d_t) * dl;
  GET_PROPERTY(particle_var::age) += d_t;

  return (GET_PROPERTY(particle_var::length) > GET_PROPERTY(particle_var::l_cp))
             ? MC::Status::Division
             : MC::Status::Idle;
  // return MC::Status::Idle;
}

void _contribution_udf([[maybe_unused]] std::size_t idx,
                       [[maybe_unused]] std::size_t position,
                       [[maybe_unused]] double weight,
                       [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
                       [[maybe_unused]] const MC::ContributionView& contribution)
{
  auto access = contribution.access();
  access(0, position) += weight * GET_PROPERTY(::particle_var::contrib_phi_s); // NOLINT
}

void _division_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                   [[maybe_unused]] std::size_t idx,
                   [[maybe_unused]] std::size_t idx2,
                   [[maybe_unused]] const MC::DynParticlesModel<float>& arr,
                   [[maybe_unused]] const MC::DynParticlesModel<float>& buffer_arr)
{
  constexpr auto local_lc = 3e-6;
  const FloatType new_current_length =
      GET_PROPERTY(particle_var::length) / static_cast<FloatType>(2.);
  GET_PROPERTY(particle_var::length) = new_current_length;
  GET_PROPERTY(particle_var::age) = 0;
  GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length) = new_current_length;
  GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::age) = 0;
  GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_cp) = local_lc;
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
  return {INDEX_FROM_ENUM(particle_var::age), INDEX_FROM_ENUM(particle_var::length)};
};

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
