#include <mc/events.hpp>


MC::EventContainer MC::EventContainer::reduce(std::span<MC::EventContainer> _data)
{
    EventContainer results;

    for (const auto& partial : _data)
    {
        for (size_t i = 0; i < n_event_type; ++i)
        {
            results.events[i] += partial.events[i];
        }
    }

    return results;
}
