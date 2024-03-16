#ifndef __MC_UNIT_HPP
#define __MC_UNIT_HPP

#include <mc/domain.hpp>

namespace MC
{
  struct MonteCarloUnit
  {
    ReactorDomain domain;

    MonteCarloUnit(const MonteCarloUnit &other) = delete;
    MonteCarloUnit(MonteCarloUnit &&other) = default;
    MonteCarloUnit() = default;
  };
} // namespace MC

#endif