#ifndef __MC_UNIT_HPP__
#define __MC_UNIT_HPP__

#include <cmt_common/macro_constructor_assignment.hpp>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/thread_private_data.hpp>


#include <variant_model.hpp>
/**
 * @namespace MC
 * @brief Namespace that contains classes and structures related to Monte Carlo
 * (MC) simulations.
 *
 * This namespace encapsulates all the components used in Monte Carlo
 * simulations, including data structures for managing simulation state and
 * results.
 */
namespace MC
{
  struct MonteCarloUnit
  {
    EventContainer events;
    ReactorDomain domain;
    AutoGenerated::ContainerVariant container;
    KPRNG rng;
    KPRNG particle_rng;
    double init_weight;
    // template <class Archive> void serialize(Archive &ar) { ar(container,domain); }
    // template <class Archive> void load(Archive &ar) { ar(container,domain); }

    SET_NON_COPYABLE(MonteCarloUnit)
    SET_DEFAULT_MOVABLE(MonteCarloUnit)
    MonteCarloUnit() = default;
    ~MonteCarloUnit() = default;
  };
} // namespace MC

#endif //__MC_UNIT_HPP__