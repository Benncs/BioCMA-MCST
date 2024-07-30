#include <mc/events.hpp>

MC::EventContainer
MC::EventContainer::reduce_local(std::span<MC::EventContainer> _data)
{
  EventContainer results;

  for (const auto &partial : _data)
  {
    for (size_t i = 0; i < n_event_type; ++i)
    {
      results.events[i] += partial.events[i];
    }
  }

  return results;
}

MC::EventContainer MC::EventContainer::reduce(std::span<size_t> _data)
{
  EventContainer results;
  for (size_t i = 0; i < _data.size(); ++i)
  {
    size_t index = i % n_event_type;
    results.events[index] += _data[i];
  }
  return results;
}