// #ifdef DECLARE_EXPORT_UDF
#ifndef __BIO_MODEL_USER_HPP__
#  define __BIO_MODEL_USER_HPP__
#include "Kokkos_Macros.hpp"
#  include <mc/traits.hpp>

namespace Models
{
  struct UdfModel
  {
  
    static std::size_t n_var;
    static constexpr std::string_view name = "udf_model";
    using uniform_weight = std::true_type; // Using type alias
    using Self = UdfModel;
    using FloatType = float;
    using SelfParticle = MC::DynParticlesModel<FloatType>;

    KOKKOS_FUNCTION static void
    init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr);

    KOKKOS_FUNCTION static double mass([[maybe_unused]] std::size_t idx,
                                              [[maybe_unused]] const SelfParticle& arr);

    KOKKOS_FUNCTION static MC::Status
    update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] FloatType d_t,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const MC::LocalConcentration& c);

    KOKKOS_FUNCTION static void
    division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
             [[maybe_unused]] std::size_t idx,
             [[maybe_unused]] std::size_t idx2,
             [[maybe_unused]] const SelfParticle& arr,
             [[maybe_unused]] const SelfParticle& buffer_arr);

    KOKKOS_FUNCTION static void
    contribution([[maybe_unused]] std::size_t idx,
                 [[maybe_unused]] std::size_t position,
                 [[maybe_unused]] double weight,
                 [[maybe_unused]] const SelfParticle& arr,
                 [[maybe_unused]] const MC::ContributionView& contributions);

    
    static std::vector<std::string_view> names();

    static std::vector<std::size_t> get_number();

    static void set_nvar();
  };
  inline std::size_t UdfModel::n_var =0;//Need to be overwritte
  static_assert(ModelType<UdfModel>, "Check Pimpl");
  static_assert(HasExportProperties<UdfModel>, "Check Pimpl");
} // namespace Models

#endif
// #endif
