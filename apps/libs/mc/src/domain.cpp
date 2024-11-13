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
  // ReactorDomain::ReactorDomain(ReactorDomain &&other) noexcept
  // {
  //   if (this != &other)
  //   {
  //     this->size = other.size;
  //     this->id = other.id;
  //     this->_total_volume = other._total_volume;
  //     this->shared_containers = other.shared_containers;
  //   }
  // }

  void ReactorDomain::setVolumes(std::span<double const> volumes_gas,
                                 std::span<double const> volumes_liq)
  {
    // Lot of assert to ensure validy of flowmap during MPI broadcast

    assert(volumes_gas.size() == size);
    assert(volumes_liq.size() == size);

    // Ok because of ShareSpace
    this->_total_volume = 0;
    for (size_t i_c = 0; i_c < volumes_gas.size(); ++i_c)
    {
      assert(volumes_liq[i_c] >= 0);
      assert(volumes_gas[i_c] >= 0);
      shared_containers(i_c).volume_liq = volumes_liq[i_c];
      shared_containers(i_c).volume_gas = volumes_gas[i_c];
      this->_total_volume += volumes_liq[i_c];
    }
  }

  ReactorDomain::ReactorDomain() : shared_containers(Kokkos::view_alloc("domain_containers"))
  {
  }

  ReactorDomain::ReactorDomain(std::span<double> volumes,
                               const CmaRead::Neighbors::Neighbors_const_view_t& _neighbors)
      : size(volumes.size()), shared_containers(Kokkos::view_alloc("domain_containers")),
        neighbors(_neighbors)
  {

    // Volume data is located on the host, creating a first unmanaged view
    Kokkos::View<double*, Kokkos::HostSpace> tmp_volume_host(volumes.data(), volumes.size());

    // Copy to data into new view that can be accessed within a kernel
    auto volume_compute =
        Kokkos::create_mirror_view_and_copy(Kokkos::SharedSpace(), tmp_volume_host);

    Kokkos::resize(shared_containers, volumes.size());

    // Temporary view for initialisation
    Kokkos::View<double, ComputeSpace> _tmp_tot("domain_tmp_total_volume", 1);
    auto tmp_shared_containers = shared_containers;
    Kokkos::parallel_for(
        "init_domain", volumes.size(), KOKKOS_LAMBDA(const int i) {
          auto& local_container = tmp_shared_containers(i);
          // Make a container with initial information about domain
          local_container = ContainerState();
          local_container.id = i;
          local_container.n_cells = 0;
          local_container.volume_liq = volume_compute(i);
          // Don't need gas volume right now, will be set during the simulation
          _tmp_tot() += local_container.volume_liq;
        });
    Kokkos::fence();

    this->_total_volume = _tmp_tot(); // copy computed volume
  }

  ReactorDomain& ReactorDomain::operator=(ReactorDomain&& other) noexcept
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

  std::vector<uint64_t> ReactorDomain::getRepartition() const
  {
    // Repartition is used for postprocessing outside kernel so we can use STL
    std::vector<uint64_t> dist(shared_containers.extent(0));

    // shared_containers is located in a sharespace so deep copy is not
    // necessary but in case of migrating shared_container to computeSpace,
    // nothing to change here
    auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), shared_containers);

    for (size_t i = 0LU; i < dist.size(); ++i)
    {
      dist[i] = static_cast<uint64_t>(host_view(i).n_cells);
    }
    return dist;
  }

  void
  ReactorDomain::in_place_reduce(std::span<const size_t> data, size_t original_size, size_t n_rank)
  {
    // OK because of sharedspace
    // Kokkos::resize(reduced.shared_containers, original_size);
    // Kokkos::View<ContainerState *, Kokkos::SharedSpace> tmp("tmp_reduce",original_size);

    if (data.size() != original_size * n_rank)
    {
      throw std::runtime_error("Cannot reduce different reactor type");
    }
    auto shared_container = this->shared_containers;
    for (size_t i_rank = 1; i_rank < n_rank; ++i_rank)
    {
      for (size_t i_c = 0; i_c < original_size; ++i_c)
      {
        shared_container(i_c).n_cells += data[i_c + i_rank * original_size];
      }
    }
  }

  ReactorDomain
  ReactorDomain::reduce(std::span<const size_t> data, size_t original_size, size_t n_rank)
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
        reduced.shared_containers(i_c).n_cells += data[i_c + i_rank * original_size];
      }
    }
    reduced.size = original_size;
    return reduced;
  }

} // namespace MC