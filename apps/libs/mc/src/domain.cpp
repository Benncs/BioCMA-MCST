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
                               const HostNeighsView& _neighbors)
      : _total_volume(std::reduce(volumes.begin(), volumes.end(), 0.)),
        size(volumes.size()),
        k_neighbor(Kokkos::view_alloc(Kokkos::WithoutInitializing, "neighbors"))

  {

    setLiquidNeighbors(_neighbors);
  }

  void ReactorDomain::setLiquidNeighbors(const HostNeighsView& data)
  {

    Kokkos::resize(k_neighbor, data.extent(0), data.extent(1));
    Kokkos::deep_copy(k_neighbor, data);
  }

  void ReactorDomain::setLiquidNeighbors(std::span<const size_t> flat_data)
  {
    const auto n_rows = this->getNumberCompartments();
    const auto n_cols = flat_data.size() / n_rows;
    KOKKOS_ASSERT(n_rows * n_cols == flat_data.size() &&
                  flat_data.size() % n_rows == 0);
    const auto* chunk = flat_data.data();
    HostNeighsView neighbors_view(chunk, n_rows, n_cols);
    setLiquidNeighbors(neighbors_view);
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
