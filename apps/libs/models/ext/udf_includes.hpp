#ifndef __BIO__EXT_MODULE_DEF__
#define __BIO__EXT_MODULE_DEF__

#include "Kokkos_Core_fwd.hpp"
#include <utility>
#ifdef DECLARE_EXPORT_UDF
#  include <dynlib/dyn_module.hpp>
#  include <dynlib/dynlib.hpp>
#  include <memory>
#endif
#include <Kokkos_Core.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
#include <models/udf_model.hpp>

/**
  @brief Unsafe namespace to handle UDF (User-defined function) via dynamic
  library loading
 */
namespace UnsafeUDF
{
  /**
  @brief Static class to access to low-level configuration to load dynamic
  library
   */
  struct Loader
  {
    static std::size_t (*set_nvar_udf)();

    static void (*init_udf)(
        const MC::KPRNG::pool_type& random_pool,
        std::size_t idx,
        const Models::UdfModel::SelfParticle& arr,
        const Models::UdfModel::Config& config); //< init function ptr

    static MC::Status (*update_udf)(
        const MC::KPRNG::pool_type& random_pool,
        float d_t,
        std::size_t idx,
        const Models::UdfModel::SelfParticle& arr,
        const MC::LocalConcentration& c); //< update function ptr

    static MC::ContribIndexBounds (*get_bounds_udf)();

    static Models::UdfModel::Config (*get_config_udf)(
        Kokkos::DefaultHostExecutionSpace& ep, std::size_t n);

    static void (*division_udf)(const MC::KPRNG::pool_type& random_pool,
                                std::size_t idx,
                                std::size_t idx2,
                                const MC::DynParticlesModel<float>& arr,
                                const MC::DynParticlesModel<float>&
                                    buffer_arr); //< division function ptr

    static double (*mass)(
        std::size_t idx,
        const MC::DynParticlesModel<float>& arr); //< mass function ptr

    static std::vector<std::string_view> (*names)();

    static std::vector<std::size_t> (*get_number)();

    // static std::vector<std::string> (*names)(); //< names function ptr

#ifdef DECLARE_EXPORT_UDF
    /**
      @brief Load UDF from .so path
     */
    [[nodiscard]] static std::shared_ptr<DynamicLibrary>
    init_lib(std::string_view path);
#endif
  };

} // namespace UnsafeUDF

#ifdef DECLARE_EXPORT_UDF

using init_udf_ptr =
    decltype(UnsafeUDF::Loader::init_udf); //< init function ptr type
using update_udf_ptr =
    decltype(UnsafeUDF::Loader::update_udf); //< update function ptr type

using get_bounds_udf_ptr = decltype(UnsafeUDF::Loader::get_bounds_udf);

using get_config_udf_ptr = decltype(UnsafeUDF::Loader::get_config_udf);

using division_udf_ptr =
    decltype(UnsafeUDF::Loader::division_udf); //< division function ptr type
using mass_udf_ptr = decltype(UnsafeUDF::Loader::mass);
using set_nvar_udf_ptr =
    decltype(UnsafeUDF::Loader::set_nvar_udf); //< division function ptr type

using names_udf_ptr = decltype(UnsafeUDF::Loader::names);
using get_number_udf_ptr =
    decltype(UnsafeUDF::Loader::get_number); //< division function ptr type

// clang-format off
/**
  @brief Module declaration
*/
DEFINE_MODULE(
              MODULE_ITEM(init_udf)
              MODULE_ITEM(update_udf)
              MODULE_ITEM(division_udf)
              MODULE_ITEM(mass_udf)
              MODULE_ITEM(names_udf)
              MODULE_ITEM(get_number_udf)
              MODULE_ITEM(set_nvar_udf)
              MODULE_ITEM(get_bounds_udf)
              MODULE_ITEM(get_config_udf)
              )
// clang-format on
#endif

#endif
