#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP


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

namespace Simulation::KernelInline
{
  struct TagFirstPass
  {
  };
  struct TagSecondPass
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

    KOKKOS_INLINE_FUNCTION
    CycleFunctor(MC::ParticlesContainer<M> _particles,
                 MC::KPRNG::pool_type _random_pool,
                 MC::KernelConcentrationType&& _concentrations,
                 kernelContribution _contribs_scatter,
                 MC::EventContainer _event)
        : d_t(0.), particles(_particles), random_pool(_random_pool),
          concentrations(std::move(_concentrations)),
          // contribs_scatter(std::move(_contribs_scatter)),
          contribs(_contribs_scatter), events(std::move(_event))
    {
    }

    void update(double _d_t, MC::ParticlesContainer<M> _particles)
    {
      this->d_t = _d_t;
      this->particles = _particles;
    }

    void sblock(MC::ContributionView _block)
    {
      this->block = _block;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const TagFirstPass _tag,
                    const TeamMember& team_handle) const
    {
      (void)_tag;
      const int idx_start = team_handle.league_rank() * team_handle.team_size();
      const int n_virtual_position = 500;
      auto access = block.access();
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team_handle, team_handle.team_size()),
          [&](const int i)
          {
            const std::size_t idx = i + idx_start;
            if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
            {
              return;
            }
            auto s = random_pool.get_state();
            const uint32_t virtual_position = s.urand(0, n_virtual_position);
            random_pool.free_state(s);
            const double weight = particles.get_weight(idx);
            const int idx_start_consumed = 2; // TODO
            const int idx_n_consumed = 1;     // TODO
            // for (int j = 0; j < idx_n_consumed; ++j)
            // {
            //   auto mock_contrib =
            //       -1.0 * weight * particles.model(idx, idx_start_consumed +
            //       j);
            //   access(virtual_position, j) += mock_contrib;
            // }
            Kokkos::parallel_for(
                Kokkos::ThreadVectorRange(team_handle, idx_n_consumed),
                [&](const int j)
                {
                  const double mock_contrib =
                      weight * particles.model(idx, idx_start_consumed + j);
                  access(virtual_position, j) += mock_contrib;
                });
          });
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TagSecondPass _tag,
                                           const TeamMember& team_handle,
                                           value_type& reduce_val) const
    {
      (void)_tag;
      (void)reduce_val.dead_total; // Counter not used currently because there
                                   // is no cell mortality
      GET_INDEX(particles.n_particles());
      if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
      {
        return;
      }

      auto local_c =
          Kokkos::subview(concentrations, Kokkos::ALL, particles.position(idx));

      if (M::update(random_pool, d_t, idx, particles.model, local_c) ==
          MC::Status::Division)
      {
        if (!particles.handle_division(random_pool, idx))
        {
          reduce_val.waiting_allocation_particle += 1;
          Kokkos::printf("[KERNEL] Division Overflow\r\n");
          events.wrap_incr<MC::EventType::Overflow>();
        }
        events.wrap_incr<MC::EventType::NewParticle>();
      };

      // particles.template get_contributions<ComputeSpace, TeamMember>(
      //     random_pool, idx, contribs_scatter, team_handle);
    }

    M::FloatType d_t;
    MC::ParticlesContainer<M> particles;
    MC::KPRNG::pool_type random_pool;
    MC::KernelConcentrationType concentrations;
    // MC::ContributionView contribs_scatter;
    kernelContribution contribs;
    MC::KernelConcentrationType limitation_factor;
    MC::EventContainer events;
    MC::ContributionView block;
  };

  template <ModelType M> struct CycleFunctorFinalReduce
  {
  public:
    CycleFunctorFinalReduce(Kokkos::View<double**, ComputeSpace> _tmp_block,
                            MC::ContributionView _contribs)
        : tmp_block(_tmp_block), contribs(_contribs)
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(const int i)
    {
      auto access = contribs.access();
      access(0, 0) += tmp_block(i, 0);
    }

    Kokkos::View<double**, ComputeSpace> tmp_block;
    MC::ContributionView contribs;
  };

} // namespace Simulation::KernelInline
#endif
