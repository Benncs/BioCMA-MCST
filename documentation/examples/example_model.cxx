#ifndef __EXAMPLE_MODEL_HPP__
#define __EXAMPLE_MODEL_HPP__

#include <common/traits.hpp>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/uptake.hpp>
#include <models/utils.hpp>
#include <string_view>

struct ModelExample
{
  // Define properties here, fill enum witn name
  enum class particle_var : int
  {
    mass = 0,
    a,
    __COUNT__ // This member is recommended to automatically get the number of
              // defined property and there get the right size
  };

  // n_var is needed, if using __COUNT__ this line does not have to be modified
  static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::__COUNT__);
  static constexpr std::string_view name =
      "example"; // Name is optional refers to model name

  using uniform_weight =
      std::true_type; //  Define this line if models works if unique weight, if
                      //  each MCparticle has
  //  specific weight comment this line
  using Self = ModelExample; // Mandatory Type alias
  using FloatType = float;   // Mandatory Type alias
  using SelfParticle =
      MC::ParticlesModel<Self::n_var, Self::FloatType>; // Mandatory Type alias

  KOKKOS_INLINE_FUNCTION static void
  init(const MC::KPRNG::pool_type& random_pool,
       std::size_t idx,
       const SelfParticle& arr);

  KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx,
                                            const SelfParticle& arr);

  KOKKOS_INLINE_FUNCTION static MC::Status
  update(const MC::KPRNG::pool_type& random_pool,
         FloatType d_t,
         std::size_t idx,
         const SelfParticle& arr,
         const MC::LocalConcentration& c);

  KOKKOS_INLINE_FUNCTION static void
  division(const MC::KPRNG::pool_type& random_pool,
           std::size_t idx,
           std::size_t idx2,
           const SelfParticle& arr,
           const SelfParticle& buffer_arr);

  KOKKOS_INLINE_FUNCTION static void
  contribution(std::size_t idx,
               std::size_t position,
               double weight,
               const SelfParticle& arr,
               const MC::ContributionView& contributions);
};

CHECK_MODEL(ModelExample)

#endif
