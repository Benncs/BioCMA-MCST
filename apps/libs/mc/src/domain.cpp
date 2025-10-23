#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <cassert>
#include <cstddef>
#include <mc/domain.hpp>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace MC
{

  void ReactorDomain::setVolumes(std::span<double const> volumes_liq)
  {

    assert(volumes_liq.size() == size);

    this->_total_volume =
        std::reduce(volumes_liq.begin(), volumes_liq.end(), 0.);

    this->_total_volume = 0;
    for (auto&& i_c : volumes_liq)
    {
      assert(i_c >= 0);
      this->_total_volume += i_c;
    }
  }

  ReactorDomain::ReactorDomain()
      : k_neighbor(Kokkos::view_alloc(Kokkos::WithoutInitializing, "neighbors"))
  {
  }

  ReactorDomain::ReactorDomain(std::span<double> volumes,
                               const NeighborsView<HostSpace>& _neighbors)
      : _total_volume(std::reduce(volumes.begin(), volumes.end(), 0.)),
        size(volumes.size()),
        k_neighbor(Kokkos::view_alloc(Kokkos::WithoutInitializing, "neighbors"))

  {

    setLiquidNeighbors(_neighbors);
  }

  void ReactorDomain::setLiquidNeighbors(const NeighborsView<HostSpace>& data)
  {

    Kokkos::resize(k_neighbor, data.extent(0), data.extent(1));
    Kokkos::deep_copy(k_neighbor, data);
  }

  ReactorDomain& ReactorDomain::operator=(ReactorDomain&& other) noexcept
  {
    if (this != &other)
    {
      this->id = other.id;
      this->size = other.size;
      this->_total_volume = other._total_volume;
    }
    return *this;
  }

  [[deprecated("Not needed anymore")]] void ReactorDomain::in_place_reduce(
      std::span<const size_t> data, size_t original_size, size_t n_rank)
  {
    // OK because of sharedspace
    // Kokkos::resize(reduced.shared_containers, original_size);
    // Kokkos::View<ContainerState *, Kokkos::SharedSpace>
    // tmp("tmp_reduce",original_size);
    if (data.size() != original_size * n_rank)
    {
      throw std::runtime_error("Cannot reduce different reactor type");
    }
    // auto shared_container = this->shared_containers;
    for (std::size_t i_rank = 1; i_rank < n_rank; ++i_rank)
    {
      for (std::size_t i_c = 0; i_c < original_size; ++i_c)
      {
        // shared_container(i_c).n_cells += data[i_c + i_rank * original_size];
      }
    }
  }

  ReactorDomain ReactorDomain::reduce(std::span<const size_t> data,
                                      size_t original_size,
                                      size_t n_rank)
  {
    ReactorDomain reduced;
    // Kokkos::resize(reduced.shared_containers, original_size);

    if (data.size() != original_size * n_rank)
    {
      throw std::runtime_error("Cannot reduce different reactor type");
    }

    for (size_t i_rank = 0; i_rank < n_rank; ++i_rank)
    {
      for (size_t i_c = 0; i_c < original_size; ++i_c)
      {
        // reduced.shared_containers(i_c).n_cells += data[i_c + i_rank *
        // original_size];
      }
    }
    reduced.size = original_size;
    return reduced;
  }

} // namespace MC