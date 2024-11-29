#include <Kokkos_Core.hpp>
#include <cassert>
#include <mc/events.hpp>
#include <span>

namespace MC
{
  EventContainer EventContainer::reduce(std::span<size_t> _data)
  {
    EventContainer results;
    assert(_data.size() % number_event_type == 0);
    for (size_t i = 0; i < _data.size(); ++i)
    {
      size_t index = i % number_event_type;
      results._events(index) += _data[i];
    }
    return results;
  }


} // namespace MC
