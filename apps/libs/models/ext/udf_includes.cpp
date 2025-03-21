#include "mc/alias.hpp"
#include "models/udf_model.hpp"
#include <mc/traits.hpp>
#ifdef DECLARE_EXPORT_UDF
#  include <dynlib/dynlib.hpp>
#  include <udf_includes.hpp>

namespace UnsafeUDF
{
  std::size_t (*Loader::set_nvar_udf)() = nullptr;

  void (*Loader::init_udf)(const MC::KPRNG::pool_type& random_pool,
                           std::size_t idx,
                           const Models::UdfModel::SelfParticle& arr) =
      nullptr; //< init function ptr

  MC::Status (*Loader::update_udf)(const MC::KPRNG::pool_type& random_pool,
                                   float d_t,
                                   std::size_t idx,
                                   const Models::UdfModel::SelfParticle& arr,
                                   const MC::LocalConcentration& c) =
      nullptr; //< update function ptr

  void (*Loader::contribution_udf)(std::size_t idx,
                                   std::size_t position,
                                   double weight,
                                   const Models::UdfModel::SelfParticle& arr,
                                   const MC::ContributionView& contributions) =
      nullptr; //< contribution function ptr

  void (*Loader::division_udf)(const MC::KPRNG::pool_type& random_pool,
                               std::size_t idx,
                               std::size_t idx2,
                               const Models::UdfModel::SelfParticle& arr,
                               const Models::UdfModel::SelfParticle& buffer_arr) =
      nullptr; //< division function ptr

  double (*Loader::mass)(std::size_t idx,
                         const Models::UdfModel::SelfParticle& arr) = nullptr; //< mass function ptr

  std::vector<std::string_view> (*Loader::names)() = nullptr;

  std::vector<std::size_t> (*Loader::get_number)() = nullptr;

  [[nodiscard]] std::shared_ptr<DynamicLibrary> Loader::init_lib(std::string_view path)
  {
    auto _handle = DynamicLibrary::getLib(path);
    auto _mod = DynamicLibrary::getModule<Module>(_handle);

    Loader::names = _mod._names_udf_m;
    Loader::get_number = _mod._get_number_udf_m;
    Loader::init_udf = _mod._init_udf_m;
    Loader::update_udf = _mod._update_udf_m;
    Loader::contribution_udf = _mod._contribution_udf_m;
    Loader::division_udf = _mod._division_udf_m;
    Loader::mass = _mod._mass_udf_m;
    Loader::set_nvar_udf = _mod._set_nvar_udf_m;
    Models::UdfModel::set_nvar();
    return _handle;
  }
} // namespace UnsafeUDF

#endif