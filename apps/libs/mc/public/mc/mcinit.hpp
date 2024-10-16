#ifndef __MC_INIT_HPP__
#define __MC_INIT_HPP__

#include <Kokkos_Atomic.hpp>
#include <Kokkos_DynamicView.hpp>
#include <cassert>
#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <mc/particles/extra_process.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>

namespace MC
{

  namespace
  {

    /**
     * @brief helper function to get initial particle weight
     */
    inline double get_initial_weight(double scale_factor,
                                     double particle_concentration,
                                     double total_volume,
                                     size_t n_particles)
    {
      // Scale factor is a fine tunning adjustement in case of misprediction of
      // particle weight particle_concentration is expected to be the real cell
      // concentration in g/L (kg/m3) Total volume is expected to be in m3 As a
      // result we can calculate the mass carried by each MC particle
      return scale_factor * particle_concentration * total_volume /
             static_cast<double>(n_particles);
    }

    /**
     */
    template <ParticleModel Model>
    void impl_init(std::unique_ptr<MonteCarloUnit> &unit,
                   double weight,
                   size_t particle_per_process)
    {

      auto container = ParticlesContainer<Model>(particle_per_process);

      auto &compartments = unit->domain.data();
      auto &list = container.get_compute();

      auto rng = unit->rng;
      auto particle_rng = list.rng_instance;

      constexpr double allocation_factor = 2.5;
      list.set_allocation_factor(allocation_factor);
      container.get_extra().extra_process.set_allocation_factor(
          allocation_factor);
      container.get_extra() =
          Results<ComputeSpace, Model>(particle_per_process);
      const auto n_compartments = unit->domain.getNumberCompartments();

      Kokkos::parallel_for(
          "mc_init",
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
              0, particle_per_process),
          KOKKOS_LAMBDA(const int i) {
            auto p = Particle<Model>(weight);
            p.properties.weight = weight;
            const uint64_t location = rng.uniform_u(0, n_compartments);
            p.properties.current_container = location;
            Kokkos::atomic_increment(&compartments(location).n_cells);
            p.init(particle_rng);
            list.set(i, std::move(p));
          });

      Kokkos::fence();

      unit->container = container;
    }

  } // namespace

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
  template <ParticleModel Model>
  std::unique_ptr<MonteCarloUnit>
  init(const ExecInfo &info,
       size_t n_particles,
       std::span<double> volumes,
       CmaRead::Neighbors::Neighbors_const_view_t neighbors,
       double x0)
  {
    auto unit = std::make_unique<MonteCarloUnit>();

    unit->domain = ReactorDomain(volumes, neighbors);

    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }

    double weight =
        get_initial_weight(1., x0, unit->domain.getTotalVolume(), n_particles);

    assert(weight > 0);
    unit->init_weight = weight;
    impl_init<Model>(unit, weight, particle_per_process);
    return unit;
  }

} // namespace MC

#endif //__MC_INIT_HPP__