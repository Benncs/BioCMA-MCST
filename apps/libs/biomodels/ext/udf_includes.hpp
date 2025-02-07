#ifndef __BIO__EXT_MODULE_DEF__
#define __BIO__EXT_MODULE_DEF__

#ifdef DECLARE_EXPORT_UDF
#  include <dynlib/dyn_module.hpp>
#  include <dynlib/dynlib.hpp>
#  include <memory>
#endif
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng.hpp>
namespace Models
{
  struct User;
  struct UserImpl;

} // namespace Models

namespace MC
{
  class ParticleDataHolder;
} // namespace MC

namespace UnsafeUDF
{
  struct Loader
  {
    static void (*make_udf)(Models::User&);
    static void (*init_udf)(Models::UserImpl&, MC::ParticleDataHolder&,MC::KPRNG);
    static void (*update_udf)(Models::UserImpl&,
                              double d_t,
                              MC::ParticleDataHolder& p,
                              const LocalConcentrationView& concentration);
    static void (*delete_udf)(Models::UserImpl**);

    static void (*contribution_udf)(Models::UserImpl&,
                                    MC::ParticleDataHolder& p,
                                    const ContributionView& contrib);

    static Models::UserImpl* (*division_udf)(Models::UserImpl&,
                                             MC::ParticleDataHolder& p,
                                             MC::KPRNG);

    static double (*mass)(Models::UserImpl&);

    static void (*fill_properties)(Models::UserImpl&, SubViewtype);

    static std::vector<std::string> (*names)();

    static std::size_t(*get_number)();

#ifdef DECLARE_EXPORT_UDF
    [[nodiscard]] static std::shared_ptr<DynamicLibrary> init_lib(std::string_view path);
#endif
  };

} // namespace UnsafeUDF

#ifdef DECLARE_EXPORT_UDF

using init_udf_ptr = decltype(UnsafeUDF::Loader::init_udf);
using make_udf_ptr = decltype(UnsafeUDF::Loader::make_udf);
using delete_udf_ptr = decltype(UnsafeUDF::Loader::delete_udf);
using update_udf_ptr = decltype(UnsafeUDF::Loader::update_udf);
using contribution_udf_ptr = decltype(UnsafeUDF::Loader::contribution_udf);
using division_udf_ptr = decltype(UnsafeUDF::Loader::division_udf);
using mass_udf_ptr = decltype(UnsafeUDF::Loader::mass);
using fill_prop_udf_ptr = decltype(UnsafeUDF::Loader::fill_properties);
using names_udf_ptr = decltype(UnsafeUDF::Loader::names);
using get_number_udf_ptr = decltype(UnsafeUDF::Loader::get_number);


// clang-format off
DEFINE_MODULE(
              MODULE_ITEM(init_udf) 
              MODULE_ITEM(make_udf) 
              MODULE_ITEM(delete_udf)
              MODULE_ITEM(update_udf)
              MODULE_ITEM(contribution_udf) 
              MODULE_ITEM(division_udf)
              MODULE_ITEM(mass_udf) 
              MODULE_ITEM(fill_prop_udf)
              MODULE_ITEM(names_udf)
              MODULE_ITEM(get_number_udf)
              )
// clang-format on
#endif

#endif