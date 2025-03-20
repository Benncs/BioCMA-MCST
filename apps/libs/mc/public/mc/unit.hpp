#ifndef __MC_UNIT_HPP__
#define __MC_UNIT_HPP__

#include "variant_model.hpp"
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
#include <models/two_meta.hpp>
#include <variant>

struct TagDetector
{
  KOKKOS_FUNCTION void operator()(const Kokkos::TeamPolicy<ComputeSpace>::member_type& team_handle,
                                  int& dead_count) const
  {
    (void)team_handle;
    (void)dead_count;
  }
  TagDetector() = default;
};

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

  template <typename FunctorType>
  Kokkos::TeamPolicy<ComputeSpace>
  get_policty(FunctorType& f, std::size_t range, bool reduce = false)
  {
    (void)f;
    (void)reduce;
    Kokkos::TeamPolicy<ComputeSpace> _policy;
    // int recommended_team_size = (reduce)
    //                                 ? _policy.team_size_recommended(f,
    //                                 Kokkos::ParallelReduceTag()) :
    //                                 _policy.team_size_recommended(f, Kokkos::ParallelForTag());

    int recommended_team_size =
        _policy.team_size_recommended(TagDetector(), Kokkos::ParallelReduceTag());
    int league_size = (static_cast<int>(range) + recommended_team_size - 1) / recommended_team_size;

    _policy = Kokkos::TeamPolicy<ComputeSpace>(league_size, recommended_team_size);

    //   std::cout << "Policy(" << league_size << "," << recommended_team_size << ")" << std::endl;
    return _policy;
  }

  /**
   * @brief General-purpose Monte Carlo unit to carry out simulations.
   *
   * The MonteCarloUnit struct encapsulates the components necessary for
   * running Monte Carlo simulations.
   * The MonteCarloUnit is non-copyable but movable, ensuring that only one
   * instance manages its resources at any time.
   */
  struct MonteCarloUnit
  {
    EventContainer events; ///< Container to manage and store simulation events
    ReactorDomain domain;  ///< Represents the domain within which the simulation occurs

    // AutoGenerated::ContainerVariant
    //     container; ///< Variant container holding various simulation entities

    AutoGenerated::ContainerVariant container;

    KPRNG rng;            ///< Random number generator used for Monte Carlo methods
    double init_weight{}; ///< Initial weight or factor used in the simulation

    MonteCarloUnit(const MonteCarloUnit&) = delete; ///< Prevent copying of the MonteCarloUnit
    MonteCarloUnit&
    operator=(const MonteCarloUnit&) = delete; ///< Prevent copying of the MonteCarloUnit
    MonteCarloUnit&
    operator=(MonteCarloUnit&&) noexcept = default;      ///< Allow moving of the MonteCarloUn
    MonteCarloUnit(MonteCarloUnit&&) noexcept = default; ///< Allow moving of the MonteCarloUn

    MonteCarloUnit() = default;
    ~MonteCarloUnit() = default;

    [[nodiscard]] uint64_t n_particle() const;

    [[nodiscard]] std::vector<uint64_t> getRepartition() const;

    template <class Archive> void serialize(Archive& ar)
    {
      ar(init_weight, events, domain, container);
      std::visit([&ar](auto& container) { ar(container); }, container);
    }
  };

} // namespace MC

#endif //__MC_UNIT_HPP__