#ifndef __MC_UNIT_HPP
#define __MC_UNIT_HPP

#include <mc/particles/particles_container.hpp>
#include "mc/prng.hpp"
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/thread_private_data.hpp>

namespace MC
{
  
  struct MonteCarloUnit
  {
    ReactorDomain domain;
    ParticlesContainer container;
    std::vector<EventContainer> ts_events;
    PRNG rand;
    std::vector<ThreadPrivateData> extras; //TODO move to unit 
    void merge(size_t i_thread){
      container.merge(extras[i_thread]);
    }
    

    MonteCarloUnit(const MonteCarloUnit &other) = delete;
    MonteCarloUnit(MonteCarloUnit &&other) = default;
    MonteCarloUnit &operator=(MonteCarloUnit &&other) = default;
    MonteCarloUnit &operator=(const MonteCarloUnit &other) = delete;
    MonteCarloUnit() = default;
    ~MonteCarloUnit() = default;
  };
} // namespace MC

#endif