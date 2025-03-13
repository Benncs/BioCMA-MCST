#ifndef __CORE_POST_PROCESS_PUBLIC_HPP__
#define __CORE_POST_PROCESS_PUBLIC_HPP__
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


//   template <typename ListType, typename ExecutionSpace>
//   void inner(std::size_t n_p,
//              ListType list,
//              ParticlePropertyViewType<ComputeSpace>& particle_values,
//              ParticlePropertyViewType<ComputeSpace>& spatial_values)
//   {
//     Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>
//         scatter_spatial_values(spatial_values);
//     Kokkos::parallel_for(
//         "kernel_get_properties",
//         Kokkos::RangePolicy<ExecutionSpace>(0, n_p),
//         KOKKOS_LAMBDA(const int i_particle) {
//           auto& p = list(i_particle);
//           auto subview = Kokkos::subview(particle_values, Kokkos::ALL, i_particle);
//           p.fill_properties(subview, scatter_spatial_values);
//         });
//     Kokkos::fence();
//     Kokkos::Experimental::contribute(spatial_values, scatter_spatial_values);
//   }

//   template <typename ListType> BonceBuffer get_properties(const std::size_t n_p ,ListType& list,const std::size_t n_compartment)
//     {
//       BonceBuffer properties;
//       // const std::size_t n_p = list.size(); // USE list size not Kokkos View size. ParticleList
//       //                                      // allocates more particles than needed

//       using ViewParticleType = ListType::value_type;
//       using Model = ViewParticleType::Model;
//       properties.vnames = Model::names();
//       properties.vnames.emplace_back("hydraulic_time");
//       properties.vnames.emplace_back("age");
//       const std::size_t n_properties = Model::get_number() + 2; //+2 for hydraulic time and age
//       // TODO: Find out if we cna use compile time size or not
//       // For the moment double** is used for simplicyt with particle model but algortihm remains
//       // exactly the same if we know the size of view.
//       // Kokkos::View<double* [n_properties], ComputeSpace> particle_values(
//       //     "device_property_values", names.size(), n_p);

//       ParticlePropertyViewType<ComputeSpace> particle_values("property_values", n_properties, n_p);
//       ParticlePropertyViewType<ComputeSpace> spatial_values(
//           "property_spatial", n_properties, n_compartment);

//       inner<ListType, Kokkos::DefaultExecutionSpace>(n_p, list, particle_values, spatial_values);

//       properties.particle_values =
//           Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particle_values);
//       properties.spatial_values =
//           Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), spatial_values);
//       return properties;
//     }



  
} // namespace PostProcessing
#endif