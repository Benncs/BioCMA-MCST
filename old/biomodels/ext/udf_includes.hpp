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
} //namespace MC

namespace UnsafeUDF
{
  // extern void (*make_udf)(Models::User &);
  // extern void (*init_udf)(Models::User::Impl &, MC::ParticleDataHolder&);
  // extern void (*delete_udf)(Models::User::Impl **);

  struct Loader
  {
    static void (*make_udf)(Models::User&);
    static void (*init_udf)(Models::UserImpl&, MC::ParticleDataHolder&);
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
DEFINE_MODULE(MODULE_ITEM(init_udf) MODULE_ITEM(make_udf) MODULE_ITEM(delete_udf)
                  MODULE_ITEM(update_udf) MODULE_ITEM(contribution_udf) MODULE_ITEM(division_udf))
#endif

#endif