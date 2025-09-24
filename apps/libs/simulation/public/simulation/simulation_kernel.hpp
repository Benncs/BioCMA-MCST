#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP
#include "Kokkos_Assert.hpp"
#include "Kokkos_Printf.hpp"
#include "simulation/move_kernel.hpp"
#include <Kokkos_Core.hpp>
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

    // KOKKOS_INLINE_FUNCTION
    // CycleFunctor(M::FloatType dt,
    // MC::ParticlesContainer<M> _particles,
    // MC::KPRNG::pool_type _random_pool,
    // MC::KernelConcentrationType&& _concentrations,
    // MC::ContributionView _contribs_scatter,
    // MC::EventContainer _event)
    //: d_t(dt), particles(_particles), random_pool(_random_pool),
    // concentrations(std::move(_concentrations)),
    // contribs_scatter(std::move(_contribs_scatter)),
    // events(std::move(_event))
    //{
    //}
    //

    KOKKOS_INLINE_FUNCTION
    CycleFunctor(MC::ParticlesContainer<M> _particles,
                 MC::KPRNG::pool_type _random_pool,
                 MC::KernelConcentrationType&& _concentrations,
                 MC::ContributionView _contribs_scatter,
                 MC::EventContainer _event)
        : d_t(0.), particles(_particles), random_pool(_random_pool),
          concentrations(std::move(_concentrations)),
          contribs_scatter(std::move(_contribs_scatter)),
          events(std::move(_event))
    {
    }

    void update(double _d_t, MC::ParticlesContainer<M> _particles)
    {
      this->d_t = _d_t;
      this->particles = _particles;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team_handle,
                                           value_type& reduce_val) const
    {
      (void)reduce_val.dead_total; // Counter not used currently because there
                                   // is no cell mortality
      GET_INDEX(particles.n_particles());
      if (particles.status(idx) != MC::Status::Idle) [[unlikely]]
      {
        // Kokkos::printf("Skip %ld", idx);
        return;
      }

      if (M::update(random_pool,
                    d_t,
                    idx,
                    particles.model,
                    Kokkos::subview(concentrations,
                                    Kokkos::ALL,
                                    particles.position(idx))) ==
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
      particles.get_contributions(idx, contribs_scatter);
    }

    M::FloatType d_t;
    MC::ParticlesContainer<M> particles;
    MC::KPRNG::pool_type random_pool;
    MC::KernelConcentrationType concentrations;
    MC::ContributionView contribs_scatter;
    MC::EventContainer events;
  };

  template <typename Space, ModelType Model> struct Functors
  {
    using FModel = Model;
    using cycle_reducer_view_type =
        KernelInline::CycleReducer<Space>::result_view_type;
    using move_reducer_view_type = Kokkos::View<std::size_t, Space>;

    typedef CycleFunctor<Model> cycle_kernel_type;
    typedef MoveFunctor move_kernel_type;
    // using cycle_kernel_type = CycleFunctor<Model>;

    cycle_reducer_view_type cycle_reducer;
    move_reducer_view_type move_reducer;

    cycle_kernel_type cycle_kernel;
    move_kernel_type move_kernel;

    explicit Functors(cycle_kernel_type&& ck, move_kernel_type&& mk)
        : cycle_reducer("cycle_reducer"), move_reducer("move_reducer"),
          cycle_kernel(std::move(ck)), move_kernel(std::move(mk))
    {
    }
  };

} // namespace Simulation::KernelInline
#endif
