// #ifdef DECLARE_EXPORT_UDF
// #include "udf_includes.hpp"
#include "ext/udf_includes.hpp"
#include <mc/alias.hpp>
#include <models/udf_model.hpp>

namespace Models
{

  KOKKOS_FUNCTION void UdfModel::init(const MC::KPRNG::pool_type& random_pool,
                                  std::size_t idx,
                                  const UdfModel::SelfParticle& arr)
  {
    UnsafeUDF::Loader::init_udf(random_pool,idx,arr);
  }

  KOKKOS_FUNCTION double UdfModel::mass(std::size_t idx, const UdfModel::SelfParticle& arr)
  {
    return UnsafeUDF::Loader::mass(idx,arr);
  }

  KOKKOS_FUNCTION MC::Status UdfModel::update(const MC::KPRNG::pool_type& random_pool,
                                          FloatType d_t,
                                          std::size_t idx,
                                          const SelfParticle& arr,
                                          const MC::LocalConcentration& c)
  {
    return UnsafeUDF::Loader::update_udf(random_pool,d_t,idx,arr,c);
  }

  KOKKOS_FUNCTION void UdfModel::division(const MC::KPRNG::pool_type& random_pool,
                                      std::size_t idx,
                                      std::size_t idx2,
                                      const SelfParticle& arr,
                                      const SelfParticle& buffer_arr)
  {
    UnsafeUDF::Loader::division_udf(random_pool,idx,idx2,arr,buffer_arr);
  }

  KOKKOS_FUNCTION void UdfModel::contribution(std::size_t idx,
                                          std::size_t position,
                                          double weight,
                                          const SelfParticle& arr,
                                          const MC::ContributionView& contributions)
  {
    UnsafeUDF::Loader::contribution_udf(idx,position,weight,arr,contributions);
  }

  void UdfModel::set_nvar()
  {
    UdfModel::n_var = UnsafeUDF::Loader::set_nvar_udf();
  }

} // namespace Models

// #endif