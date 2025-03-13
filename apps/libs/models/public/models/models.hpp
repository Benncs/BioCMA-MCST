#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include <mc/traits.hpp>

namespace Models
{

  struct SimpleModel
  {
    enum particle_var
    {
      mass = 0,
    };
    static constexpr std::size_t n_var = 1;
    static constexpr std::string_view name = "simple";
    using Self = SimpleModel;
    using FloatType = float;
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    KOKKOS_INLINE_FUNCTION static void init(std::size_t idx, const SelfParticle& arr)
    {
      arr(idx, particle_var::mass) = 1.;
    }

    KOKKOS_INLINE_FUNCTION static MC::Status
    update(FloatType d_t, std::size_t idx, const SelfParticle& arr, const MC::LocalConcentration& c)
    {
      (void)c;
      arr(idx, particle_var::mass) += d_t * 1;
      // return (idx % 3 == 0) ? MC::Status::Division : MC::Status::Idle;
      return MC::Status::Idle;
    }

    KOKKOS_INLINE_FUNCTION static void division(std::size_t idx,
                                                std::size_t idx2,
                                                const SelfParticle& arr,
                                                const SelfParticle& buffer_arr)
    {
      // for(size_t i=0;i<n_var;++i)
      // {
      //   buffer_arr(idx,i)=arr(idx,i); //basic copy
      // }
      const FloatType new_current_mass = arr(idx, particle_var::mass) / static_cast<FloatType>(2.);
      buffer_arr(idx2, particle_var::mass) = new_current_mass;
      arr(idx, particle_var::mass) = new_current_mass;
    }
  };

  CHECK_MODEL(SimpleModel)

} // namespace Models

#endif