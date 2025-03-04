#ifndef __SIMULATION_MC_KERNEL_UPDATE_HPP
#define __SIMULATION_MC_KERNEL_UPDATE_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <common/kokkos_vector.hpp>
#include <mc/container_state.hpp>
#include <mc/events.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/mcparticles.hpp>
#include <simulation/alias.hpp>

static constexpr bool const_number_simulation = false;
namespace Simulation::KernelInline
{
  using GeneratorPool = Kokkos::Random_XorShift64_Pool<>; // TODO remove
  template <typename ListType> class Kernel
  {

  public:
    KOKKOS_FUNCTION void operator()(std::size_t i_particle) const;
    Kernel(double _d_t,
           ListType _list,
           ListType _extra,
           Kokkos::View<MC::ContainerState*, Kokkos::SharedSpace> _local_compartments,
           MC::EventContainer _events,
           ContributionView _biomass_contribution,
           Kokkos::View<size_t, Kokkos::SharedSpace> _waiting_allocation_particle,
           GeneratorPool pool)
        : d_t(_d_t), list(_list), extra(std::move(_extra)),
          local_compartments(std::move(_local_compartments)), events(std::move(_events)),
          biomass_contribution(std::move(_biomass_contribution)),
          waiting_allocation_particle(std::move(_waiting_allocation_particle)), rng_pool(pool)
    {
    }

  private:
    double d_t{};
    ListType list;
    ListType extra;
    Kokkos::View<MC::ContainerState*, Kokkos::SharedSpace> local_compartments;
    GeneratorPool rng_pool;
    MC::EventContainer events;
    ContributionView biomass_contribution;
    Kokkos::View<size_t, Kokkos::SharedSpace> waiting_allocation_particle;

    KOKKOS_INLINE_FUNCTION void handle_division(auto& particle) const
    {

      if (particle.properties.status == MC::CellStatus::CYTOKINESIS)
      {
        events.template incr<MC::EventType::NewParticle>();

        auto new_particle = particle.division(list.rng_instance);
        KOKKOS_ASSERT(new_particle.properties.status != MC::CellStatus::CYTOKINESIS);

        if constexpr (!const_number_simulation)
        {

          if (!extra.emplace(std::move(new_particle)))
          {
            Kokkos::printf("SPAWNING OVERFLOW\r\n");
            Kokkos::atomic_increment(&waiting_allocation_particle()); // TODO clean this
          }
        }
        else
        {
          particle.properties.weight *= 2;
          Kokkos::atomic_increment(
              &local_compartments(particle.properties.current_container).n_cells);
        }

        KOKKOS_ASSERT(particle.properties.status != MC::CellStatus::CYTOKINESIS);
      }
    }
  };

  /****************
  ** IMPL KERNEL **
  *****************/

  template <typename ListType>
  KOKKOS_FUNCTION void Kernel<ListType>::operator()(const std::size_t i_particle) const
  {

    auto& particle = list._owned_data[i_particle];
    auto& status = particle.properties.status;
    auto& properties = particle.properties;

    // // // TODO Test attribute
    if (status == MC::CellStatus::DEAD) [[unlikely]]
    {
      // This warp divergence is not a priori a problem because it is a
      // cold branch
      return;
    }

    particle.update(d_t, local_compartments(properties.current_container).concentrations, rng_pool);

    particle.contribution(biomass_contribution);

    handle_division(particle);
  }

} // namespace Simulation::KernelInline

#endif //!__SIMULATION_MC_KERNEL_UPDATE_HPP