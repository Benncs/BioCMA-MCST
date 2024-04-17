#ifndef __SIMULATION_VEC_MAPFLOWS_HPP__
#define __SIMULATION_VEC_MAPFLOWS_HPP__
#include <vector>

#include <simulation/matflows.hpp>

namespace Simulation
{
  
  struct BasicCacheMatflows
  {
    std::vector<MatFlow> data;

    explicit BasicCacheMatflows(size_t n) noexcept
    {
      data.resize(n);
    }
  };
} // namespace Simulation

#endif //__SIMULATION_VEC_MAPFLOWS_HPP__
