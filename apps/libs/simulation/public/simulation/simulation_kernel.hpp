#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP
#include "mc/particles/particles_container.hpp"
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
                                        const MC::pool_type& _random_pool,
                                        MC::KernelConcentrationType _concentrations)
        : d_t(dt), particles(_particles), random_pool(_random_pool),
          concentrations(std::move(_concentrations))
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team_handle,
                                           int& waiting_allocation_particle,
                                           int& dead_total) const
    {
      (void)dead_total; // Counter not used currently because there is no cell mortality
      GET_INDEX(particles.n_particles());
      if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
      {
        return;
      }

      if (M::update(d_t,
                    idx,
                    particles.model,
                    Kokkos::subview(concentrations, particles.position(idx), Kokkos::ALL)) ==
          MC::Status::Division)
      {
        if (!particles.handle_division(idx))
        {
          waiting_allocation_particle += 1;
          Kokkos::printf("Division Overflow\r\n");
          assert(false && "Division Overflow Not implemented");
        }
      };
    }

    M::FloatType d_t;
    MC::ParticlesContainer<M> particles;
    MC::pool_type random_pool;
    MC::KernelConcentrationType concentrations;
  };

} // namespace Simulation::KernelInline
#endif
