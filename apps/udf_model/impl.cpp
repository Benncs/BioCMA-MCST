#include "mc/macros.hpp"
#include <cstdlib>
#include <mc/traits.hpp>
#include <udf_includes.hpp>
using namespace Models;

void __attribute__((constructor)) on_load()
{
  printf("UDF loaded\n"); // NOLINT
}

enum class particle_var : int
{
  mass = 0,
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
  GET_PROPERTY(particle_var::mass) = 1.e-15 * gen.frand();
  random_pool.free_state(gen);
};

MC::Status _update_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                       [[maybe_unused]] float d_t,
                       [[maybe_unused]] std::size_t idx,
                       [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
                       [[maybe_unused]] const MC::LocalConcentration& c)
{

  // NOP
  return MC::Status::Idle;
}

void _contribution_udf([[maybe_unused]] std::size_t idx,
                       [[maybe_unused]] std::size_t position,
                       [[maybe_unused]] double weight,
                       [[maybe_unused]] const Models::UdfModel::SelfParticle& arr,
                       [[maybe_unused]] const MC::ContributionView& contribution)
{
}

void _division_udf([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                   [[maybe_unused]] std::size_t idx,
                   [[maybe_unused]] std::size_t idx2,
                   [[maybe_unused]] const MC::DynParticlesModel<float>& arr,
                   [[maybe_unused]] const MC::DynParticlesModel<float>& buffer_arr)
{
}

double mass(std::size_t idx, const Models::UdfModel::SelfParticle& arr)
{
  return GET_PROPERTY(particle_var::mass);
}

// clang-format off
EXPORT_MODULE(
    module, &_init_udf, 
    &_update_udf, 
    &_contribution_udf, 
    &_division_udf, 
    &mass, 
    &_set_nvar);
//clang-format on