#include <common/kokkos_vector.hpp>
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <cstring>
#include <mc/unit.hpp>
#include <utility>
#include <variant>

template <typename ListType, typename ScatterType> struct NcellFunctor
{
  ScatterType n_cells;
  ListType list;

  KOKKOS_FUNCTION NcellFunctor(ListType _list, ScatterType _n_cells)
      : n_cells(std::move(_n_cells)), list(std::move(_list))
  {
  }
  KOKKOS_FUNCTION void operator()(const int i_particle) const
  {
    auto access_ = n_cells.access();
    access_(list._owned_data(i_particle).properties.current_container) += 1;
  };
};

namespace MC
{
  [[nodiscard]] uint64_t MonteCarloUnit::n_particle()const
  {
    return std::visit([](auto&& c) -> std::size_t { return c.n_particle(); }, container);
  }

  [[nodiscard]] std::vector<uint64_t> MonteCarloUnit::getRepartition()const
  {
    if (domain.getNumberCompartments() == 1)
    {
      return {n_particle()};
    }
    // Compute on-deman ncells as long as we don´t need it during simulation but only while
    // exporting. If ncells has to be known more often (id during iteration), n_cell view can be
    // allocated in mc_unit and incremented during kernel using atomic.
    Kokkos::View<uint64_t*, ComputeSpace> n_cells("n_cell", this->domain.getNumberCompartments());
    Kokkos::Experimental::ScatterView<uint64_t*> sn_cells(n_cells);
    std::visit(
        [&sn_cells](auto&& _container)
        {
          auto& list = _container.get_compute();
          Kokkos::parallel_for("get_repartition", list.size(), NcellFunctor(list, sn_cells));
        },
        container);
    Kokkos::Experimental::contribute(n_cells, sn_cells);
    auto host = Kokkos::create_mirror_view_and_copy(HostSpace(), n_cells);
    std::vector<uint64_t> dist(n_cells.extent(0));
    std::memcpy(dist.data(), host.data(), host.size() * sizeof(uint64_t));

    return dist;
  }
} // namespace MC