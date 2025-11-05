#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP

#include "mc/alias.hpp"
#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Printf.hpp>
#include <Kokkos_Random.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles_container.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
#include <simulation/alias.hpp>
#include <simulation/probability_leaving.hpp>
#include <simulation/probe.hpp>
#include <utility>

#define CHECK_STATUS_OR_RETURN(__idx__)                                        \
  if (particles.status(__idx__) != MC::Status::Idle) [[unlikely]]              \
  {                                                                            \
    return;                                                                    \
  }

namespace Simulation::KernelInline
{
  struct TagContribution
  {
  };
  struct TagCycle
  {
  };

  struct CycleReduceType
  {
    std::size_t waiting_allocation_particle;
    std::size_t dead_total;
  };

  template <class Space> class CycleReducer
  {
  public:
    // Required for Concept
    using reducer = CycleReducer;
    using value_type = CycleReduceType;
    using result_view_type = Kokkos::View<value_type, Space>;

    KOKKOS_INLINE_FUNCTION
    void join(value_type& dest, const value_type& src) const
    {
      dest.dead_total += src.dead_total;
      dest.waiting_allocation_particle += src.waiting_allocation_particle;
    }

    KOKKOS_INLINE_FUNCTION
    value_type& reference() const
    {
      return *value.data();
    }

    KOKKOS_INLINE_FUNCTION
    result_view_type view() const
    {
      return value;
    }

    KOKKOS_INLINE_FUNCTION
    bool references_scalar() const
    {
      return references_scalar_v;
    }

    // Optional
    KOKKOS_INLINE_FUNCTION
    void init(value_type& val) const
    {
      val.dead_total = 0;
      val.waiting_allocation_particle = 0;
    }

    // KOKKOS_INLINE_FUNCTION
    // void final(value_type& val) const
    // {
    //   // NOP
    // }

    // Part of Build-In reducers for Kokkos
    KOKKOS_INLINE_FUNCTION
    CycleReducer(value_type& value_) : value(&value_), references_scalar_v(true)
    {
    }

    KOKKOS_INLINE_FUNCTION
    CycleReducer(const result_view_type& value_)
        : value(value_), references_scalar_v(false)
    {
    }

  private:
    result_view_type value;
    bool references_scalar_v;
  };
  template <ModelType M> struct CycleFunctor
  {
    using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
    using TeamMember = TeamPolicy::member_type;
    using value_type = CycleReduceType;

    CycleFunctor() = default;

    KOKKOS_INLINE_FUNCTION
    CycleFunctor(MC::ParticlesContainer<M> _particles,
                 MC::KPRNG::pool_type _random_pool,
                 MC::KernelConcentrationType&& _concentrations,
                 MC::ContributionView _contribs_scatter,
                 MC::EventContainer _event)
        : d_t(0.), particles(std::move(_particles)), random_pool(_random_pool),
          concentrations(std::move(_concentrations)),
          contribs_scatter(std::move(_contribs_scatter)),
          events(std::move(_event))
    {
    }

    void update(double _d_t, MC::ParticlesContainer<M> _particles)
    {
      this->d_t = _d_t;
      this->particles = std::move(_particles);
    }

    // KOKKOS_INLINE_FUNCTION
    // void operator()(const TagFirstPass _tag,
    //                 const TeamMember& team_handle) const
    // {
    //   (void)_tag;
    //   GET_INDEX(particles.n_particles());
    //   if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
    //   {
    //     return;
    //   }
    //   particles.get_contributions(idx, contribs_scatter);
    // }
    //
    KOKKOS_INLINE_FUNCTION
    void operator()(const TagContribution _tag, const std::size_t idx) const
    {
      (void)_tag;
      CHECK_STATUS_OR_RETURN(idx);

      particles.get_contributions(idx, contribs_scatter);
    }

    // KOKKOS_INLINE_FUNCTION void operator()(const TagCycle _tag,
    //                                        const TeamMember& team_handle,
    //                                        value_type& reduce_val) const
    // {

    //   (void)_tag;
    //   (void)reduce_val.dead_total; // Counter not used currently because
    //                                // there  is no cell mortality
    //   // GET_INDEX(particles.n_particles());

    //   std ::size_t start_idx =
    //       (team_handle.league_rank() * team_handle.team_size());
    //   std ::size_t idx = start_idx + team_handle.team_rank();
    //   // if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
    //   // {
    //   //   return;
    //   // }
    //   const std::size_t n_species = concentrations.extent(0);
    //   Kokkos::View<double**, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    //       shmem_view(
    //           team_handle.team_shmem(), team_handle.team_size(), n_species);

    //   Kokkos::parallel_for(
    //       Kokkos::TeamThreadRange(team_handle, 0, team_handle.team_size()),
    //       [&](int i)
    //       {
    //         std::size_t particle_idx = start_idx + i;
    //         const auto position = particles.position(particle_idx);
    //         for (int s = 0; s < n_species; ++s)
    //         {
    //           shmem_view(i, s) = concentrations(s, position);
    //         }
    //       });
    //   // team_handle.team_barrier();
    //   CHECK_STATUS_OR_RETURN(idx);
    //   if (idx > particles.n_particles())
    //   {
    //     return;
    //   }
    //   particles.ages(idx, 1) += d_t;

    //   auto local_c =
    //       Kokkos::subview(shmem_view, team_handle.team_rank(), Kokkos::ALL);

    //   if (M::update(random_pool, d_t, idx, particles.model, local_c) ==
    //       MC::Status::Division)
    //   {
    //     if (!particles.handle_division(random_pool, idx))
    //     {
    //       reduce_val.waiting_allocation_particle += 1;
    //       events.wrap_incr<MC::EventType::Overflow>();
    //       Kokkos::printf("[KERNEL] Division Overflow\r\n");
    //     }
    //     events.wrap_incr<MC::EventType::NewParticle>();
    //   };
    // }

    KOKKOS_INLINE_FUNCTION void
    operator()(TagCycle _tag,
               const Kokkos::TeamPolicy<ComputeSpace>::member_type& team,
               value_type& reduce_val) const
    {
      (void)_tag;
      const auto team_size = team.team_size();

      const auto work_stride = team_size * team.league_size();
      const auto work_start = team.league_rank() * team_size + team.team_rank();

      for (auto idx = work_start; idx < particles.n_particles();
           idx += work_stride)
      {
        CHECK_STATUS_OR_RETURN(idx);
        exec_per_particle(idx, reduce_val);
      }
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TagCycle _tag,
                                           const std::size_t idx,
                                           value_type& reduce_val) const
    {

      (void)_tag;
      (void)reduce_val.dead_total;
      exec_per_particle(idx, reduce_val);
    }

    KOKKOS_INLINE_FUNCTION void exec_per_particle(const std::size_t idx,
                                                  value_type& reduce_val) const
    {
      particles.ages(idx, 1) += d_t;
      auto local_c =
          Kokkos::subview(concentrations, Kokkos::ALL, particles.position(idx));

      if (M::update(random_pool, d_t, idx, particles.model, local_c) ==
          MC::Status::Division)
      {
        if (!particles.handle_division(random_pool, idx))
        {
          reduce_val.waiting_allocation_particle += 1;
          events.wrap_incr<MC::EventType::Overflow>();
          Kokkos::printf("[KERNEL] Division Overflow\r\n");
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
    // kernelContribution contribs;
    MC::KernelConcentrationType limitation_factor;
    MC::EventContainer events;
  };

} // namespace Simulation::KernelInline
#endif
