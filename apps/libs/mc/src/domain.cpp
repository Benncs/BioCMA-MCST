#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <cassert>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <mc/container_state.hpp>
#include <mc/domain.hpp>
#include <stdexcept>
#include <utility>

namespace
{
  template <typename ViewType, typename VolumeViewType> class FunctorInit
  {
  public:
    FunctorInit(ViewType shared_containers, VolumeViewType _volumes)
        : tmp_shared_containers(std::move(shared_containers)), volumes(std::move(_volumes))
    {
    }

    ViewType tmp_shared_containers;
    VolumeViewType volumes;

    KOKKOS_FUNCTION void operator()(std::size_t i_particle,double& total_volume) const
    {
      auto& local_container = tmp_shared_containers(i_particle);
      // Make a container with initial information about domain
      local_container = MC::ContainerState();
      local_container.id = i_particle;
      local_container.n_cells = 0;
      local_container.volume_liq = volumes(i_particle);
      // Don't need gas volume right now, will be set during the simulation
      total_volume+=local_container.volume_liq;
    }
  };
} // namespace

namespace MC
{

  void ReactorDomain::setVolumes(std::span<double const> volumes_gas,
                                 std::span<double const> volumes_liq)
  {

    assert(volumes_gas.size() == size);
    assert(volumes_liq.size() == size);

    // Ok because of ShareSpace
    this->_total_volume = 0;
    for (size_t i_c = 0; i_c < volumes_gas.size(); ++i_c)
    {
      assert(volumes_liq[i_c] >= 0);
      assert(volumes_gas[i_c] >= 0);
      shared_containers[i_c].volume_liq = volumes_liq[i_c];
      shared_containers[i_c].volume_gas = volumes_gas[i_c];
      
      this->_total_volume += volumes_liq[i_c];
    }
  }

  ReactorDomain::ReactorDomain()
      : shared_containers(Kokkos::view_alloc("domain_containers")),
        k_neighbor(Kokkos::view_alloc(Kokkos::WithoutInitializing, "neighbors"))
  {
  }

  ReactorDomain::ReactorDomain(std::span<double> volumes,
                               const NeighborsView<HostSpace>& _neighbors)
      : size(volumes.size()), shared_containers(Kokkos::view_alloc("domain_containers")),
        k_neighbor(Kokkos::view_alloc(Kokkos::WithoutInitializing, "neighbors"))

  {
    setLiquidNeighbors(_neighbors);
    // Volume data is located on the host, creating a first unmanaged view
    Kokkos::View<double*, Kokkos::HostSpace> tmp_volume_host(volumes.data(), volumes.size());

    // Copy to data into new view that can be accessed within a kernel
    auto volume_compute =
        Kokkos::create_mirror_view_and_copy(Kokkos::SharedSpace(), tmp_volume_host);

    Kokkos::resize(shared_containers, volumes.size());

    // Temporary view for initialisation
    auto tmp_shared_containers = shared_containers;
    double local_total_volume=0;
    Kokkos::parallel_reduce("init_domain",
                         Kokkos::RangePolicy<ComputeSpace>(0, volumes.size()),
                         FunctorInit(shared_containers, volume_compute),Kokkos::Sum<double>(local_total_volume));

    Kokkos::fence();
    this->_total_volume = local_total_volume;
  }

  void ReactorDomain::setLiquidNeighbors(const NeighborsView<HostSpace>& data)
  {

    Kokkos::resize(k_neighbor, data.extent(0), data.extent(1));
    Kokkos::deep_copy(k_neighbor,data);

  }

  ReactorDomain& ReactorDomain::operator=(ReactorDomain&& other) noexcept
  {
    if (this != &other)
    {
      this->id = other.id;
      this->size = other.size;
      // this->neighbors = other.neighbors;
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

  [[deprecated("Not needed anymore")]] void
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
    for (std::size_t i_rank = 1; i_rank < n_rank; ++i_rank)
    {
      for (std::size_t i_c = 0; i_c < original_size; ++i_c)
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