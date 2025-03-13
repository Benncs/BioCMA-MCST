#ifndef __CORE_POST_PROCESS_PUBLIC_HPP__
#define __CORE_POST_PROCESS_PUBLIC_HPP__
#include "mc/particles/particles_container.hpp"
#include "mc/traits.hpp"
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <mc/particles/particle_model.hpp>

template <typename MemorySpace>
using ParticlePropertyViewType = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

using SubViewtype =
    Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>, decltype(Kokkos::ALL), std::size_t>;
namespace PostProcessing
{
  struct BonceBuffer
  {
    ParticlePropertyViewType<HostSpace> particle_values;
    ParticlePropertyViewType<HostSpace> spatial_values;
    std::vector<std::string> vnames;
  };

  template <HasExportProperties Model, typename ExecutionSpace>
  void inner(std::size_t n_p,
             const MC::ParticlePositions& position,
             const typename Model::SelfParticle& model,
             ParticlePropertyViewType<ComputeSpace>& particle_values,
             ParticlePropertyViewType<ComputeSpace>& spatial_values)
  {
    Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>
        scatter_spatial_values(spatial_values);
    Kokkos::parallel_for(
        "kernel_spatial_get_properties",
        Kokkos::RangePolicy<ExecutionSpace>(0, n_p),
        KOKKOS_LAMBDA(const int i_particle) {
          auto access = scatter_spatial_values.access();

          for (std::size_t i = 0; i < Model::n_var; ++i)
          {
            access(position(i_particle), i) += model(i_particle, i);
            particle_values( i,i_particle) = model(i_particle, i);
          }
        });
    Kokkos::fence();
    Kokkos::Experimental::contribute(spatial_values, scatter_spatial_values);
  }

  template <ModelType M>
  std::optional<PostProcessing::BonceBuffer>
  get_properties(const MC::ParticlesContainer<M>& container, const std::size_t n_compartment)
  {
    if constexpr (HasExportProperties<M>)
    {
      BonceBuffer properties;
      const std::size_t n_p =
          container.n_particles(); // USE list size not Kokkos View size. ParticleList
                                   // allocates more particles than needed

      properties.vnames = M::names();
      ParticlePropertyViewType<ComputeSpace> spatial_values(
          "property_spatial", M::n_var, n_compartment);
      ParticlePropertyViewType<ComputeSpace> particle_values("property_values", M::n_var, n_p);

      inner<M, Kokkos::DefaultExecutionSpace>(
          n_p, container.position, container.model, particle_values, spatial_values);

      properties.particle_values =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particle_values);
      properties.spatial_values =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), spatial_values);
      return properties;
    }
    else
    {
      return std::nullopt;
    }
  }

} // namespace PostProcessing
#endif