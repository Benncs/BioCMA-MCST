#ifndef __MC_EVENTS_HPP__
#define __MC_EVENTS_HPP__

#include <array>
#include <cstddef>
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

  template <EventType event> consteval size_t eventIndex()
  {
    return static_cast<size_t>(event);
  }

  struct EventContainer
  {
    std::array<size_t, n_event_type> events{0, 0, 0,0,0};

    template <EventType event> [[nodiscard]] constexpr size_t get() const
    {
      return events[eventIndex<event>()];
    }

    template <EventType event> [[nodiscard]] constexpr size_t &get_mut()
    {
      return events[eventIndex<event>()];
    }

    template <EventType event> constexpr void incr()
    {
      events[eventIndex<event>()]++;
    }

    static EventContainer reduce_local(std::span<EventContainer> _data);
    static EventContainer reduce(std::span<size_t> _data);
  };

}; // namespace MC

#endif //__MC_EVENTS_HPP__