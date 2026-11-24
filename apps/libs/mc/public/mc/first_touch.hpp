#ifndef __MC_FIRST_TOUCH_HPP__
#define __MC_FIRST_TOUCH_HPP__

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <biocma_cst_config.hpp>
#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <common/kokkos_getpolicy.hpp>
#include <cstddef>
#include <utility>

/** @file NUMA first touch for the particle arrays.

Linux places a physical page on the NUMA node of the thread that first
writes it : first touch friendly

Touching the live range explicitly, with the same team decomposition the move
kernel uses, puts each page on the domain that will actually read it. The tail
[n_live, n_allocated) is touched too it must be initialised, since
Kokkos::resize would have zeroed it but its partition does not matter until
particles grow into it.
*/
namespace MC
{
  namespace Impl
  {
    template <typename ViewType> struct FirstTouchFunctor
    {
      using TeamMember = typename Kokkos::TeamPolicy<ComputeSpace>::member_type;
      using value_type = typename ViewType::non_const_value_type;

      ViewType view;
      std::size_t begin;
      std::size_t end;
      std::size_t per_team;

      KOKKOS_INLINE_FUNCTION void
      operator()(const TeamMember& team) const
      {
        const std::size_t p0 = begin + (team.league_rank() * per_team);
        const std::size_t upper = (p0 + per_team > end) ? end - p0 : per_team;

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, upper),
                             [&](const std::size_t idx)
                             {
                               const std::size_t i = p0 + idx;
                               if constexpr (ViewType::rank == 1)
                               {
                                 view(i) = value_type{};
                               }
                               else
                               {
                                 for (std::size_t j = 0; j < view.extent(1);
                                      ++j)
                                 {
                                   view(i, j) = value_type{};
                                 }
                               }
                             });
      }
    };
  } // namespace Impl

  /// Write-touch [begin, end) with the move kernel's team decomposition.
  template <typename ViewType>
  void
  first_touch(const ViewType& view,
              std::size_t begin,
              std::size_t end,
              std::size_t npt)
  {
    end = std::min(end, static_cast<std::size_t>(view.extent(0)));
    if (begin >= end || npt == 0)
    {
      return;
    }

    const std::size_t league = Common::c_league_size(end - begin, npt);

    Kokkos::parallel_for(
        "first_touch",
        Kokkos::TeamPolicy<ComputeSpace>(
            static_cast<int>(league), Kokkos::AUTO(), Kokkos::AUTO()),
        Impl::FirstTouchFunctor<ViewType>{ view, begin, end, npt });
  }

  /** @brief Kokkos::resize, but with the first touch placed by us.

  Same observable result as Kokkos::resize -- contents preserved up to the
  overlap, remainder zeroed but the pages of the live range are placed by
  the threads that will read them.

  Like Kokkos::resize, this is a no-op when the leading extent already matches:
  reallocating would throw away the placement it is here to establish.

  @param n_live the range the compute kernels will actually walk, i.e. the
                used particle count, not the allocation.
  @param npt    particles per team, from KernelDispatchOptions::m_p_p_team_move.
                It must be what CycleFunctors launches the move kernel with,
                otherwise the touch and the compute walk the array with
                different partitions and the whole point is lost.
  */
  template <typename ViewType, typename... Extents>
  void
  resize_first_touch(ViewType& view,
                     std::size_t n_alloc,
                     std::size_t n_live,
                     std::size_t npt,
                     Extents... rest)
  {
    if (view.extent(0) == n_alloc)
    {
      return;
    }

    ViewType next(Kokkos::view_alloc(Kokkos::WithoutInitializing, view.label()),
                  n_alloc,
                  rest...);

    const std::size_t live = std::min(n_live, n_alloc);
    first_touch(next, 0, live, npt);       // placed with the compute partition
    first_touch(next, live, n_alloc, npt); // must exist; placement not yet
                                           // relevant

    const std::size_t n_copy
        = std::min(static_cast<std::size_t>(view.extent(0)), n_alloc);
    if (n_copy != 0)
    {
      const auto range = Kokkos::make_pair(std::size_t{ 0 }, n_copy);
      if constexpr (ViewType::rank == 1)
      {
        Kokkos::deep_copy(Kokkos::subview(next, range),
                          Kokkos::subview(view, range));
      }
      else
      {
        Kokkos::deep_copy(Kokkos::subview(next, range, Kokkos::ALL),
                          Kokkos::subview(view, range, Kokkos::ALL));
      }
    }

    view = next;
  }
} // namespace MC

#endif //__MC_FIRST_TOUCH_HPP__
