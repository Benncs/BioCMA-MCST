#include "Kokkos_Assert.hpp"
#include "common/common.hpp"
#include "mc/alias.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <cassert>
#include <mc/domain.hpp>
#include <numeric>

namespace MC
{

  void
  ReactorDomain::setVolumes(std::span<double const> volumes_liq)
  {

    assert(volumes_liq.size() == size);

    this->_total_volume
        = std::reduce(volumes_liq.begin(), volumes_liq.end(), 0.);
    Kokkos::View<const double*, HostSpace> tmp_host_volume(volumes_liq.data(),
                                                           volumes_liq.size());
    Kokkos::deep_copy(this->inner.liquid_volume, tmp_host_volume);
  }

  ReactorDomain::ReactorDomain(double total_volume, std::size_t size)
      : _total_volume(total_volume), size(size)
  {
  }

  ReactorDomain::ReactorDomain() : ReactorDomain(0, 0)

  {
  }

  ReactorDomain::ReactorDomain(std::span<double> volumes)
      : ReactorDomain(std::reduce(volumes.begin(), volumes.end(), 0.),
                      volumes.size())

  {
  }

  void
  ReactorDomain::update(std::span<const double> newliquid_volume,
                        std::span<const std::size_t> neighors_flat,
                        std::span<const double> out_flows,
                        std::span<const double> proba_flat)
  {

    const auto n_rows = this->getNumberCompartments();
    const auto n_cols = neighors_flat.size() / n_rows;
    if (proba_flat.size() != neighors_flat.size())
    {
      throw std::invalid_argument(
          "Neighbors and proba should have the same size");
    }

    this->setLiquidNeighbors(n_rows, n_cols, neighors_flat);
    this->setVolumes(newliquid_volume);

    KOKKOS_ASSERT(newliquid_volume.size()
                  == this->inner.liquid_volume.extent(0));

    KOKKOS_ASSERT(proba_flat.size() % n_rows == 0);

    DiagonalView<HostSpace, true> _diag_transition(out_flows.data(), n_rows);
    Kokkos::deep_copy(this->inner.diag_transition, _diag_transition);

    const auto* chunk_proba = proba_flat.data();
    CumulativeProbabilityView<HostSpace, true> tmp_host_proba(
        chunk_proba, n_rows, n_cols);
    Kokkos::resize(this->inner.cumulative_probability, n_rows, n_cols);
    Kokkos::deep_copy(this->inner.cumulative_probability, tmp_host_proba);
  }

  void
  ReactorDomain::setLiquidNeighbors(const std::size_t e1,
                                    const std::size_t e2,
                                    std::span<const size_t> flat_data)
  {

    KOKKOS_ASSERT(e1 * e2 == flat_data.size() && flat_data.size() % e1 == 0);

    using HostNeighsView = Kokkos::View<
        const std::size_t**,
        Kokkos::LayoutRight,
        HostSpace,
        Kokkos::MemoryTraits<Kokkos::RandomAccess | Kokkos::Unmanaged>>;

    const auto* chunk = flat_data.data();
    HostNeighsView neighbors_view(chunk, e1, e2);
    Kokkos::resize(this->inner.neighbors, e1, e2);
    Kokkos::deep_copy(this->inner.neighbors, neighbors_view);
  }

  void
  ReactorDomain::set_leaving_flow(const std::size_t i,
                                  const std::size_t i_flow,
                                  const double flow,
                                  const double volume) const
  {
    KOKKOS_ASSERT(flow > 0 && volume > 0);

    this->inner.leaving_flow(i)
        = { .index = i_flow, .flow = flow, .volume = volume };
  }

  ReactorDomain&
  ReactorDomain::operator=(ReactorDomain&& other) noexcept
  {
    if (this != &other)
    {
      this->id = other.id;
      this->size = other.size;
      this->_total_volume = other._total_volume;
    }
    return *this;
  }

  void
  ReactorDomain::init_inner(const std::size_t n_flows)
  {
    constexpr bool is_const = false;
    // inner = MoveInfo<ComputeSpace,false>(this->getNumberCompartments(),
    // n_flows);
    const auto n_compartments = this->getNumberCompartments();

    // clang-format off
    auto _inner = DomainState<ComputeSpace,is_const>{  };
    _inner.neighbors = MC::NeighborsView<ComputeSpace,is_const>("neighbors", 0, 0);
    _inner.diag_transition = MC::DiagonalView<ComputeSpace,is_const>("diag_transition", n_compartments);
    _inner.leaving_flow = MC::LeavingFlowView<is_const>("leaving_flow", n_flows);
    _inner.liquid_volume = MC::VolumeView<ComputeSpace,is_const>("liquid_volume", n_compartments);
    _inner.cumulative_probability = MC::CumulativeProbabilityView<ComputeSpace,is_const>("cumulative_proba",  0, 0);
    // clang-format on
    this->inner = _inner;
  }

} // namespace MC
