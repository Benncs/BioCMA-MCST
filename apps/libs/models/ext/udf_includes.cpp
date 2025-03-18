#include <mc/traits.hpp>
#ifdef DECLARE_EXPORT_UDF
#  include <dynlib/dynlib.hpp>
#  include <udf_includes.hpp>

namespace UnsafeUDF
{
  void (*Loader::make_udf)(Models::User&) = nullptr;
  void (*Loader::init_udf)(Models::UserImpl&, MC::ParticleDataHolder&, MC::KPRNG) = nullptr;
  void (*Loader::delete_udf)(Models::UserImpl**) = nullptr;
  void (*Loader::update_udf)(Models::UserImpl&,
                             double d_t,
                             MC::ParticleDataHolder& p,
                             const LocalConcentrationView& concentration) = nullptr;

  void (*Loader::contribution_udf)(Models::UserImpl&,
                                   MC::ParticleDataHolder& p,
                                   const ContributionView& contrib) = nullptr;

  double (*Loader::mass)(Models::UserImpl&) = nullptr;

  void (*Loader::fill_properties)(Models::UserImpl&, SubViewtype) = nullptr;

  Models::UserImpl* (*Loader::division_udf)(Models::UserImpl&,
                                            MC::ParticleDataHolder& p,
                                            MC::KPRNG) = nullptr;

  std::vector<std::string> (*Loader::names)() = nullptr;

  std::size_t (*Loader::get_number)() = nullptr;

  [[nodiscard]] std::shared_ptr<DynamicLibrary> Loader::init_lib(std::string_view path)
  {
    auto _handle = DynamicLibrary::getLib(path);
    auto _mod = DynamicLibrary::getModule<Module>(_handle);

    Loader::make_udf = _mod._make_udf_m;
    Loader::init_udf = _mod._init_udf_m;
    Loader::delete_udf = _mod._delete_udf_m;
    Loader::update_udf = _mod._update_udf_m;
    Loader::contribution_udf = _mod._contribution_udf_m;
    Loader::division_udf = _mod._division_udf_m;
    Loader::mass = _mod._mass_udf_m;
    Loader::fill_properties = _mod._fill_prop_udf_m;
    Loader::get_number = _mod._get_number_udf_m;
    Loader::names = _mod._names_udf_m;
    // Loader::names = _mod.nam
    return _handle;
  }
} // namespace UnsafeUDF

#endif