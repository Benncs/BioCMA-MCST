#include "common/kokkos_vector.hpp"
#include "mc/container_state.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <cassert>
#include <cma_read/neighbors.hpp>
#include <cstddef>
#include <mc/domain.hpp>
#include <stdexcept>

namespace MC
{
  ReactorDomain::ReactorDomain(ReactorDomain &&other) noexcept
  {
    if (this != &other)
    {
      this->size = other.size;
      this->id = other.id;
      this->_total_volume = other._total_volume;
      this->shared_containers = other.shared_containers;
    }
  }

  void ReactorDomain::setVolumes(std::span<double const> volumes_gas,
                                 std::span<double const> volumes_liq)
  {
    //Lot of assert to ensure validy of flowmap during MPI broadcast 
    
    assert(volumes_gas.size()==size);
    assert(volumes_liq.size()==size);

    // Ok because of ShareSpace
    this->_total_volume = 0;
    for (size_t i_c = 0; i_c < volumes_gas.size(); ++i_c)
    {
      assert(volumes_liq[i_c]>=0);
      assert(volumes_gas[i_c]>=0);
      shared_containers(i_c).volume_liq = volumes_liq[i_c];
      shared_containers(i_c).volume_gas = volumes_gas[i_c];
      this->_total_volume += volumes_liq[i_c];
    }
  }
  ReactorDomain::ReactorDomain(
      std::span<double> volumes,
      const CmaRead::Neighbors::Neighbors_const_view_t &_neighbors)
      : size(volumes.size()), neighbors(_neighbors)
  {


    Kokkos::View<double *, Kokkos::HostSpace> tmp_volume_host(volumes.data(),
                                                              volumes.size());

    auto volume_compute =
        Kokkos::create_mirror_view_and_copy(ComputeSpace(), tmp_volume_host);

    auto tmp_shared_containers =
        Kokkos::View<ContainerState *, Kokkos::SharedSpace>("domain_containers",
                                                            volumes.size());
    Kokkos::View<double, ComputeSpace> _tmp_tot("domain_tmp_total_volume", 1);

    Kokkos::parallel_for(
        "init_domain", volumes.size(), KOKKOS_LAMBDA(const int i) {
          auto &local_container = tmp_shared_containers(i);
          local_container = ContainerState();
          local_container.id = i;
          local_container.n_cells=0;
          local_container.volume_liq = volume_compute(i);
          _tmp_tot() += local_container.volume_liq;
        });
    Kokkos::fence();
    shared_containers = tmp_shared_containers;

    this->_total_volume = _tmp_tot();
  }

  ReactorDomain &ReactorDomain::operator=(ReactorDomain &&other) noexcept
  {
    if (this != &other)
    {
      this->id = other.id;
      this->size = other.size;
      this->neighbors = other.neighbors;
      this->_total_volume = other._total_volume;
      this->shared_containers = other.shared_containers;
    }
    return *this;
  }

  std::vector<size_t> ReactorDomain::getRepartition() const
  {

    std::vector<size_t> dist(shared_containers.extent(0));

    auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                         shared_containers);

    for (size_t i = 0LU; i < dist.size(); ++i)
    {
      dist[i] = static_cast<size_t>(host_view(i).n_cells);
    }
    return dist;
  }

  ReactorDomain ReactorDomain::reduce(std::span<const size_t> data,
                                      size_t original_size,
                                      size_t n_rank)
  {
    // OK because of sharedspace
    ReactorDomain reduced;
    Kokkos::resize(reduced.shared_containers, original_size);

    if (data.size() != original_size * n_rank)
    {
      throw std::runtime_error("Cannot reduce different reactor type");
    }

    for (size_t i_rank = 0; i_rank < n_rank; ++i_rank)
    {
      for (size_t i_c = 0; i_c < original_size; ++i_c)
      {
        reduced.shared_containers(i_c).n_cells +=
            data[i_c + i_rank * original_size];
      }
    }

    return reduced;
  }

  // ReactorDomain ReactorDomain::reduce(std::span<ReactorDomain> data)
  // {
  //   ReactorDomain reduced;
  //   size_t size = data[0].containers.size();
  //   auto& containers = reduced.containers;
  //   for(auto&& domain : data)
  //   {
  //     if(size!=domain.n_compartments())
  //     {
  //       throw std::runtime_error("Cannot reduce diffferent reactor type");
  //     }
  //     for(size_t i =0 ;i<size;++i)
  //     {
  //       containers[i].n_cells+=domain.containers[i].n_cells;
  //     }
  //   }

  //   return reduced;
  // }

} // namespace MC