#ifndef __MC_UNIT_HPP__
#define __MC_UNIT_HPP__

#include <cmt_common/macro_constructor_assignment.hpp>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/thread_private_data.hpp>


namespace MC
{

  struct MonteCarloUnit
  {
    ReactorDomain domain;
    ParticlesContainer container;
    MC::EventContainer events;
    std::vector<MC::Particles> extra_process;
    // std::vector<EventContainer> ts_events; //TODO Move to thread private 
    PRNG rand;
    // std::vector<ThreadPrivateData> extras; 
    // void merge(size_t i_thread)
    // {
    //   container.merge(extras[i_thread]);
    // }

    template <class Archive> void serialize(Archive &ar) { ar(container,domain); }
    template <class Archive> void load(Archive &ar) { ar(container,domain); }

    SET_NON_COPYABLE(MonteCarloUnit)
    SET_DEFAULT_MOVABLE(MonteCarloUnit)
    MonteCarloUnit() = default;
    ~MonteCarloUnit() = default;
  };
} // namespace MC

#endif //__MC_UNIT_HPP__