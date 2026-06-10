#ifndef __SIMULATION_TIMES_HPP__
#define __SIMULATION_TIMES_HPP__

namespace Simulation
{
  struct SimulatimeTimes
  {
    // t0 is 0 if new simulation or non zero for loaded simulation
    double starting_time{ 0. };

    // current relative simulation time
    double relative_time{ 0. };

    inline void
    advance(const double dt)
    {
      relative_time += dt;
    }

    [[nodiscard]] inline double
    relative() const
    {
      return relative_time;
    }

    [[nodiscard]] inline double
    absolute() const
    {
      return starting_time + relative_time;
    }
  };
} // namespace Simulation

#endif
