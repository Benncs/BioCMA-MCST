// #ifdef DECLARE_EXPORT_UDF
// #include "udf_includes.hpp"
#include "common/common.hpp"
#include "ext/udf_includes.hpp"
#include <Kokkos_Core_fwd.hpp>
#include <fwd/Kokkos_Fwd_OPENMP.hpp>
#include <mc/alias.hpp>
#include <models/udf_model.hpp>

namespace Models
{

  KOKKOS_FUNCTION void UdfModel::init(const MC::KPRNG::pool_type& random_pool,
                                      std::size_t idx,
                                      const UdfModel::SelfParticle& arr,
                                      const UdfModel::Config& config)
  {
    UnsafeUDF::Loader::init_udf(random_pool, idx, arr, config);
  }

  KOKKOS_FUNCTION double UdfModel::mass(std::size_t idx,
                                        const UdfModel::SelfParticle& arr)
  {
    return UnsafeUDF::Loader::mass(idx, arr);
  }

  KOKKOS_FUNCTION MC::Status
  UdfModel::update(const MC::KPRNG::pool_type& random_pool,
                   FloatType d_t,
                   std::size_t idx,
                   const SelfParticle& arr,
                   const MC::LocalConcentration& c)
  {
    return UnsafeUDF::Loader::update_udf(random_pool, d_t, idx, arr, c);
  }

  KOKKOS_FUNCTION void
  UdfModel::division(const MC::KPRNG::pool_type& random_pool,
                     std::size_t idx,
                     std::size_t idx2,
                     const SelfParticle& arr,
                     const SelfParticle& buffer_arr)
  {
    UnsafeUDF::Loader::division_udf(random_pool, idx, idx2, arr, buffer_arr);
  }

  void UdfModel::set_nvar()
  {
    UdfModel::n_var = UnsafeUDF::Loader::set_nvar_udf();
  }

  std::vector<std::string_view> UdfModel::names()
  {
    return UnsafeUDF::Loader::names();
  }

  MC::ContribIndexBounds UdfModel::get_bounds()
  {
    return UnsafeUDF::Loader::get_bounds_udf();
  }

  UdfModel::Config UdfModel::get_config(std::size_t n)
  {
    auto e = MC::HostSpace();
    return UnsafeUDF::Loader::get_config_udf(e, n);
  }

  std::vector<std::size_t> UdfModel::get_number()
  {
    return UnsafeUDF::Loader::get_number();
  }

} // namespace Models

// #endif
