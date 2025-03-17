#ifndef __MC_INIT_HPP__
#define __MC_INIT_HPP__

#include "mc/traits.hpp"
#include <cassert>
#include <common/execinfo.hpp>
#include <cstdint>
#include <mc/domain.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <utility>

namespace MC
{

  void impl_init(double& total_mass,
                 uint64_t n_particles,
                 MonteCarloUnit& unit,
                 ContainerVariant&& container);
  /**
   * @brief Helper function to initialize a MonteCarloUnit.
   *
   * Since MonteCarloUnit is not a generic type and the model type is resolved
   * at runtime, no constructors are defined to avoid carrying template
   * functions when using the unit. This function wraps the constructor
   * externally, providing a convenient way to initialize a MonteCarloUnit with
   * the appropriate model type and parameters.
   *
   * @tparam Model The particle model type, specified at compile time.
   *
   *
   * @return A unique pointer to the initialized MonteCarloUnit.
   */
  template <ModelType Model>
  std::unique_ptr<MonteCarloUnit> init(uint64_t n_particles,
                                       std::span<double> volumes,
                                       const NeighborsView<HostSpace>& neighbors,
                                       double& total_mass)
  {
    if constexpr (ConstWeightModelType<Model>)
    {
      Kokkos::printf("Const Weights\r\n");
    }
    auto unit = std::make_unique<MonteCarloUnit>();
    unit->domain = ReactorDomain(volumes, neighbors);
    auto container = ParticlesContainer<Model>(n_particles);
    impl_init(total_mass, n_particles, *unit, std::move(container));
    return unit;
  }

  void post_init_weight(std::unique_ptr<MonteCarloUnit>& unit, double x0, double total_mass);

} // namespace MC

#endif //__MC_INIT_HPP__
