#ifndef __MC_UNIT_HPP
#define __MC_UNIT_HPP

#include <mc/domain.hpp>
#include <mc/events.hpp>

namespace MC
{

  struct MonteCarloUnit
  {
    ReactorDomain domain;
    std::vector<EventContainer> ts_events;

    MonteCarloUnit(const MonteCarloUnit &other) = delete;
    MonteCarloUnit(MonteCarloUnit &&other) = default;
    MonteCarloUnit &operator=(MonteCarloUnit &&other) = default;
    MonteCarloUnit &operator=(const MonteCarloUnit &other) = delete;
    MonteCarloUnit() = default;
    ~MonteCarloUnit() = default;
  };
} // namespace MC

#endif