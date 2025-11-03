#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <cassert>
#include <mc/domain.hpp>
#include <numeric>

namespace MC
{

  void ReactorDomain::setVolumes(std::span<double const> volumes_liq)
  {

    assert(volumes_liq.size() == size);

    this->_total_volume =
        std::reduce(volumes_liq.begin(), volumes_liq.end(), 0.);

    this->_total_volume = 0;
    // for (auto&& i_c : volumes_liq)
    // {
    //   assert(i_c >= 0);
    //   this->_total_volume += i_c;
    // }
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

} // namespace MC
