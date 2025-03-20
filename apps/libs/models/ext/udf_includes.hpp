#ifndef __BIO__EXT_MODULE_DEF__
#define __BIO__EXT_MODULE_DEF__

#ifdef DECLARE_EXPORT_UDF
#  include <dynlib/dyn_module.hpp>
#  include <dynlib/dynlib.hpp>
#  include <memory>
#endif
#include <Kokkos_Core.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>

namespace Models
{
  struct User;
  struct UserImpl;

} // namespace Models

namespace MC
{
  class ParticleDataHolder;
} // namespace MC

/**
  @brief Unsafe namespace to handle UDF (User-defined function) via dynamic library loading
 */
namespace UnsafeUDF
{
  /**
  @brief Static class to access to low-level configuration to load dynamic library
   */
  struct Loader
  {
    static void (*make_udf)(Models::User&); //< make function ptr
    static void (*init_udf)(Models::UserImpl&,
                            MC::ParticleDataHolder&,
                            MC::KPRNG); //< init function ptr
    static void (*update_udf)(Models::UserImpl&,
                              double d_t,
                              MC::ParticleDataHolder& p,
                              const MC::LocalConcentration& concentration); //< update function ptr
    static void (*delete_udf)(Models::UserImpl**);                          //< delete function ptr

    static void (*contribution_udf)(
        Models::UserImpl&,
        MC::ParticleDataHolder& p,
        const MC::ContributionView& contrib); //< contribution function ptr

    static Models::UserImpl* (*division_udf)(Models::UserImpl&,
                                             MC::ParticleDataHolder& p,
                                             MC::KPRNG); //< division function ptr

    static double (*mass)(Models::UserImpl&); //< mass function ptr

    static std::vector<std::string> (*names)(); //< names function ptr

#ifdef DECLARE_EXPORT_UDF
    /**
      @brief Load UDF from .so path
     */
    [[nodiscard]] static std::shared_ptr<DynamicLibrary> init_lib(std::string_view path);
#endif
  };

} // namespace UnsafeUDF

#ifdef DECLARE_EXPORT_UDF

using init_udf_ptr = decltype(UnsafeUDF::Loader::init_udf);     //< init function ptr type
using make_udf_ptr = decltype(UnsafeUDF::Loader::make_udf);     //< make function ptr type
using delete_udf_ptr = decltype(UnsafeUDF::Loader::delete_udf); //< delete function ptr type
using update_udf_ptr = decltype(UnsafeUDF::Loader::update_udf); //< update function ptr type
using contribution_udf_ptr =
    decltype(UnsafeUDF::Loader::contribution_udf); //< contribution function ptr type
using division_udf_ptr = decltype(UnsafeUDF::Loader::division_udf); //< division function ptr type
using mass_udf_ptr = decltype(UnsafeUDF::Loader::mass);

using names_udf_ptr = decltype(UnsafeUDF::Loader::names); //< names function ptr type

// clang-format off
/**
  @brief Module declaration
*/ 
DEFINE_MODULE(
              MODULE_ITEM(init_udf) 
              MODULE_ITEM(make_udf) 
              MODULE_ITEM(delete_udf)
              MODULE_ITEM(update_udf)
              MODULE_ITEM(contribution_udf) 
              MODULE_ITEM(division_udf)
              MODULE_ITEM(mass_udf) 
      
              MODULE_ITEM(names_udf)
       
              )
// clang-format on
#endif

#endif