#ifndef __MC_EVENTS_HPP__
#define __MC_EVENTS_HPP__

#include "common/kokkos_vector.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Macros.hpp>
#include <array>
#include <common/execinfo.hpp>
#include <cstddef>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <span>

namespace MC
{
  enum class EventType : char
  {
    ChangeWeight = 0,
    NewParticle = 1,
    Death = 2,
    Move = 3,
    Exit = 4
  };

  constexpr size_t n_event_type = 5;
  constexpr size_t number_event_type = n_event_type;
  template <EventType event> consteval size_t eventIndex()
  {
    return static_cast<size_t>(event);
  }

  struct alignas(ExecInfo::cache_line_size) EventContainer
  {

    Kokkos::View<size_t[n_event_type], Kokkos::SharedHostPinnedSpace>
        _events; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    EventContainer()
    {
      _events = Kokkos::View<size_t[n_event_type],
                             Kokkos::SharedHostPinnedSpace>(
          "events"); // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
      Kokkos::deep_copy(_events, 0);
    }

    [[nodiscard]] std::span<size_t> get_span() const
    {
      return {_events.data(), number_event_type};
    }

    void clear() const
    {
      Kokkos::deep_copy(_events, 0);
    }

    void inplace_reduce(const EventContainer &other) const
    {

      for (size_t i = 0LU; i < n_event_type; ++i)
      {
        _events(i) += other._events(i);
      }
    }

    template <EventType event> [[nodiscard]] constexpr size_t get() const
    {
      return _events(eventIndex<event>());
    }

    template <EventType event> [[nodiscard]] constexpr size_t &get_mut()
    {
      return _events(eventIndex<event>());
    }

    template <EventType event> constexpr void incr() const
    {
      Kokkos::atomic_increment(&_events(eventIndex<event>()));
    }

    static EventContainer reduce_local(std::span<EventContainer> _data);

    static EventContainer reduce(std::span<size_t> _data);
  };
}; // namespace MC

#endif //__MC_EVENTS_HPP__