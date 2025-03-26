#ifndef __CORE_POST_PROCESS_PUBLIC_HPP__
#define __CORE_POST_PROCESS_PUBLIC_HPP__
#include "Kokkos_Core_fwd.hpp"
#include "common/common.hpp"
#include "mc/traits.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <mc/particles_container.hpp>

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

  template <HasExportPropertiesFull Model, typename ExecutionSpace>
  void inner(std::size_t n_p,
             const MC::ParticlePositions& position,
             const typename Model::SelfParticle& model,
             ParticlePropertyViewType<ComputeSpace>& particle_values,
             ParticlePropertyViewType<ComputeSpace>& spatial_values)
  {
    Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>
        scatter_spatial_values(spatial_values);

    Kokkos::parallel_for(
        "kernel_get_properties",
        Kokkos::RangePolicy<ExecutionSpace>(0, n_p),
        KOKKOS_LAMBDA(const int i_particle) {
          auto access = scatter_spatial_values.access();

          for (std::size_t i = 0; i < Model::n_var; ++i)
          {
            access(i, position(i_particle)) += model(i_particle, i);
            particle_values(i, i_particle) = model(i_particle, i);
          }
          auto mass = Model::mass(i_particle, model);
          access(Model::n_var, position(i_particle)) += mass;
          particle_values(Model::n_var, i_particle) = mass;
        });
    Kokkos::fence();
    Kokkos::Experimental::contribute(spatial_values, scatter_spatial_values);
  }

  template <HasExportPropertiesPartial Model, typename ExecutionSpace>
  void inner_partial(std::size_t n_p,
                     const MC::ParticlePositions& position,
                     const typename Model::SelfParticle& model,
                     ParticlePropertyViewType<ComputeSpace>& particle_values,
                     ParticlePropertyViewType<ComputeSpace>& spatial_values)
  {
    Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>
        scatter_spatial_values(spatial_values);
    static const auto indices = Model::get_number();

    static Kokkos::View<size_t*,HostSpace> host_index("host_index",indices.size());

    for(auto i =0;i< indices.size();++i)
    {
      host_index(i)=indices[i];
    }

    const Kokkos::View<const size_t*,ComputeSpace>  kindices = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_index);
   
    Kokkos::parallel_for(
        "kernel_get_properties",
        Kokkos::RangePolicy<ExecutionSpace>(0, n_p),
        KOKKOS_LAMBDA(const int i_particle) {
          auto access = scatter_spatial_values.access();
       
          for (int i = 0; i < kindices.extent(0); ++i)
          {
            auto index_export = kindices(i);
            access(i, position(i_particle)) += model(i_particle, index_export);
            particle_values(i, i_particle) = model(i_particle, index_export);
          }
          auto mass = Model::mass(i_particle, model);
          access(kindices.extent(0), position(i_particle)) += mass;
          particle_values(kindices.extent(0), i_particle) = mass;
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
      auto ar = M::names();
      properties.vnames = std::vector<std::string>(ar.begin(), ar.end());
      properties.vnames.emplace_back("mass");

      size_t n_var = M::n_var + 1;
      if constexpr (HasExportPropertiesPartial<M>)
      {
        n_var = ar.size() + 1;
      }
      ParticlePropertyViewType<ComputeSpace> spatial_values(
          "property_spatial", n_var, n_compartment);
      ParticlePropertyViewType<ComputeSpace> particle_values("property_values", n_var, n_p);

      if constexpr (HasExportPropertiesPartial<M>)
      {
        
        inner_partial<M, Kokkos::DefaultExecutionSpace>(
            n_p, container.position, container.model, particle_values, spatial_values);
      }
      else
      {

        inner<M, Kokkos::DefaultExecutionSpace>(
            n_p, container.position, container.model, particle_values, spatial_values);
      }

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