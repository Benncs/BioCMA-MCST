#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP
#include "Kokkos_Assert.hpp"
#include "Kokkos_Printf.hpp"
#include <mc/particles_container.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
#include <simulation/alias.hpp>
#include <simulation/probability_leaving.hpp>
#include <simulation/probe.hpp>
#include <utility>

namespace Simulation::KernelInline
{

  template <ModelType M> struct CycleFunctor
  {
    using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
    using TeamMember = TeamPolicy::member_type;

    KOKKOS_INLINE_FUNCTION CycleFunctor(M::FloatType dt,
                                        MC::ParticlesContainer<M> _particles,
                                        MC::KPRNG::pool_type _random_pool,
                                        MC::KernelConcentrationType&& _concentrations,
                                        MC::ContributionView _contribs_scatter,
                                        MC::EventContainer _event)
        : d_t(dt), particles(_particles), random_pool(_random_pool),
          concentrations(std::move(_concentrations)),
          contribs_scatter(std::move(_contribs_scatter)), events(std::move(_event))
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team_handle,
                                           std::size_t& waiting_allocation_particle,
                                           [[maybe_unused]] std::size_t& dead_total) const
    {
      (void)dead_total; // Counter not used currently because there is no cell mortality
      GET_INDEX(particles.n_particles());
      if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
      {
        Kokkos::printf("Skip %ld", idx);
        return;
      }

      if (M::update(random_pool,
                    d_t,
                    idx,
                    particles.model,
                    Kokkos::subview(concentrations, particles.position(idx), Kokkos::ALL)) ==
          MC::Status::Division)
      {
        if (!particles.handle_division(random_pool, idx))
        {
          waiting_allocation_particle += 1;
          Kokkos::printf("Division Overflow\r\n");
          events.wrap_incr<MC::EventType::Overflow>();

          // KOKKOS_ASSERT(false && "Division Overflow Not implemented");
        }
        events.wrap_incr<MC::EventType::NewParticle>();
      };

      particles.get_contributions(idx, contribs_scatter);
    }

    M::FloatType d_t;
    MC::ParticlesContainer<M> particles;
    MC::KPRNG::pool_type random_pool;
    MC::KernelConcentrationType concentrations;
    MC::ContributionView contribs_scatter;
    MC::EventContainer events;
  };

} // namespace Simulation::KernelInline
#endif
