#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP

#include "mc/particles/data_holder.hpp"
#include "simulation/probe.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <common/kokkos_vector.hpp>
#include <mc/container_state.hpp>
#include <mc/events.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>

static constexpr bool const_number_simulation = false;

#include <simulation/alias.hpp>

namespace
{

  // KOKKOS_INLINE_FUNCTION float _ln(float x)
  // {
  //   unsigned int bx = *reinterpret_cast<unsigned int *>(&x);
  //   const unsigned int ex = bx >> 23;
  //   const signed int t = static_cast<signed int>(ex) - static_cast<signed int>(127);
  //   unsigned int s = (t < 0) ? (-t) : t;
  //   bx = 1065353216 | (bx & 8388607);
  //   x = *reinterpret_cast<float *>(&bx);
  //   return -1.49278 + (2.11263 + (-0.729104 + 0.10969 * x) * x) * x + 0.6931471806 * t;
  // }


  KOKKOS_INLINE_FUNCTION bool probability_leaving(double random_number,
                                                  double volume,
                                                  double flow,
                                                  double dt)
  {
    // Use of Kokkos log to be sure that log is well defined
    // return (dt * flow / volume) > (-Kokkos::log(1 - random_number));

    // We can use 1-random_number if we have enough particle

    // return (dt * flow / volume) > (-ln(random_number));

    return (dt * flow / volume) > (-Kokkos::log(random_number));

    // const double x = volume/flow;
    // const double ax = dt/(volume/flow);
    // return random_number<ax*(1-ax/2.);
  }
} // namespace

namespace Simulation::KernelInline
{
  template <typename ListType, typename ResultViewType>
  class Kernel
  {

  public:
    KOKKOS_FUNCTION void operator()(std::size_t i_particle) const;
    Kernel(double d_t,
           ListType list,
           ResultViewType rview,
           Kokkos::View<size_t, Kokkos::SharedSpace> _internal_counter_dead,
           Kokkos::View<MC::ContainerState *, Kokkos::SharedSpace>
               _local_compartments,
           NeighborsViewCompute _neighbors,
           MC::KPRNG _local_rng,
           DiagonalViewCompute _diag_transition,
           CumulativeProbabilityViewCompute _cumulative_probability,
           MC::EventContainer _events,
           ContributionView _biomass_contribution,
           LeavingFlowType _local_leaving_flow,
           LeavingFlowIndexType _local_index_leaving_flow,
           Probes prob)
        : d_t(d_t), list(list), rview(rview),
          internal_counter_dead(std::move(_internal_counter_dead)),
          local_compartments(std::move(_local_compartments)),
          neighbors(std::move(_neighbors)), local_rng(std::move(_local_rng)),
          diag_transition(std::move(_diag_transition)),
          cumulative_probability(std::move(_cumulative_probability)),
          events(std::move(_events)),
          biomass_contribution(std::move(_biomass_contribution)),
          local_leaving_flow(std::move(_local_leaving_flow)),
          local_index_leaving_flow(std::move(_local_index_leaving_flow)),
          probe_at_t(std::move(prob))
    {
    }

  private:
    double d_t{};
    ListType list;
    ResultViewType rview;
    Kokkos::View<size_t, Kokkos::SharedSpace> internal_counter_dead;
    Kokkos::View<MC::ContainerState *, Kokkos::SharedSpace> local_compartments;
    NeighborsViewCompute neighbors;
    MC::KPRNG local_rng;
    DiagonalViewCompute diag_transition;
    CumulativeProbabilityViewCompute cumulative_probability;
    MC::EventContainer events;
    ContributionView biomass_contribution;
    LeavingFlowType local_leaving_flow;
    LeavingFlowIndexType local_index_leaving_flow;
    Probes probe_at_t;

    KOKKOS_INLINE_FUNCTION void handle_division(auto &particle) const;

    KOKKOS_INLINE_FUNCTION void handle_move(MC::ParticleDataHolder &properties) const;

    KOKKOS_INLINE_FUNCTION void handle_exit(auto &status,
                                            MC::ParticleDataHolder &properties) const;

    KOKKOS_INLINE_FUNCTION size_t __find_next_compartment(
        std::size_t i_compartment,
        double random_number,
        const Kokkos::View<const size_t *, Kokkos::LayoutStride, ComputeSpace>
            &i_neighbor) const;
/*
    KOKKOS_INLINE_FUNCTION size_t _find_next_compartment_2(
        const std::size_t i_compartment,
        const double random_number,
        const NeighborsViewCompute& n) const;i*/
  };

  /****************
  ** IMPL KERNEL **
  *****************/

  template <typename ListType, typename ResultViewType>
  KOKKOS_FUNCTION void Kernel<ListType, ResultViewType>::operator()(
      const std::size_t i_particle) const
  {

    auto &particle = list._owned_data[i_particle];
    auto &status = particle.properties.status;
    auto& properties = particle.properties;



    // // // TODO Test attribute
    if (status == MC::CellStatus::DEAD) [[unlikely]]
    {
      // This warp divergence is not a priori a problem because it is a
      // cold branch
      return;
    }

    handle_move(properties);
    handle_exit(status, properties);

    // // Useless to try to reducing branhcing here because we're going to
    // // early return so serializing will occur
    if (status == MC::CellStatus::OUT)
    {
      events.template incr<MC::EventType::Exit>();
      status = MC::CellStatus::DEAD;
      Kokkos::atomic_increment(&internal_counter_dead());
      return;
    }

    particle.update(d_t,
                    local_compartments(properties.current_container)
                        .concentrations,
                    list.rng_instance);

    particle.contribution(biomass_contribution);

    handle_division(particle);
  }

  template <typename ListType, typename ResultViewType>
  KOKKOS_INLINE_FUNCTION void
  Kernel<ListType, ResultViewType>::handle_division(auto &particle) const
  {
    if (particle.properties.status == MC::CellStatus::CYTOKINESIS)
    {
      events.template incr<MC::EventType::NewParticle>();

      const auto new_particle = particle.division(list.rng_instance);
      if constexpr (!const_number_simulation)
      {
        auto *const np = rview().extra_process.spawn();
        if (np != nullptr)
        {
          *np = new_particle;
        }
        else [[unlikely]] // TODO Test attribute
        {
          Kokkos::printf("SPAWNING OVERFLOW\r\n");
          Kokkos::atomic_increment(&rview().waiting_allocation_particle);
        }
      }
      else
      {
        particle.properties.weight *= 2;
        Kokkos::atomic_increment(
            &local_compartments(particle.properties.current_container).n_cells);
      }

      KOKKOS_ASSERT(particle.properties.status != MC::CellStatus::CYTOKINESIS);
      KOKKOS_ASSERT(new_particle.properties.status !=
                    MC::CellStatus::CYTOKINESIS);
    }
  }
  template <typename ListType, typename ResultViewType>
  KOKKOS_INLINE_FUNCTION void
  Kernel<ListType, ResultViewType>::handle_move(MC::ParticleDataHolder &properties) const
  {
    const size_t i_compartment = properties.current_container;
    auto &current_container = local_compartments(i_compartment);

    const auto i_neighbor =
        Kokkos::subview(neighbors, i_compartment, Kokkos::ALL);

    // Need 2 random numbers, use index 0 and 1 to acess to
    // const auto random = local_rng.template double_uniform<2>();

    auto generator = local_rng.random_pool.get_state();
    const double rng1 = generator.drand(0., 1.);
    local_rng.random_pool.free_state(generator);

    if (probability_leaving(rng1,
                            current_container.volume_liq,
                            diag_transition(i_compartment),
                            d_t))
    {
        auto generator = local_rng.random_pool.get_state();
        const double rng2 = generator.drand(0., 1.);
        local_rng.random_pool.free_state(generator);
      // Find the next compartment based on the random number and
      // cumulative
      // probabilities
      // const size_t next =
      //     _find_next_compartment_2(i_compartment, rng2, neighbors);
      const size_t next =
          __find_next_compartment(i_compartment, rng2, i_neighbor);

      KOKKOS_ASSERT(Kokkos::atomic_load(&current_container.n_cells) > 0); // ??
      Kokkos::atomic_decrement(&current_container.n_cells);
      Kokkos::atomic_increment(&local_compartments(next).n_cells);
      properties.current_container = next;
      events.template incr<MC::EventType::Move>();
    }
  };

  template <typename ListType, typename ResultViewType>
  KOKKOS_INLINE_FUNCTION void
  Kernel<ListType, ResultViewType>::handle_exit(auto &status,
                                                MC::ParticleDataHolder &properties) const
  {
    for (size_t i = 0LU; i < local_index_leaving_flow.size(); ++i)
    {
      const double random_number = local_rng.double_uniform();
      const auto &index = local_index_leaving_flow(i);
      const auto &flow = local_leaving_flow(i);

      if (properties.current_container != index ||
          status != MC::CellStatus::IDLE)
      {
        return;
      }

      auto &current_container = local_compartments(index);
      if (probability_leaving(
              random_number, current_container.volume_liq, flow, d_t))
      {

        if constexpr (AutoGenerated::FlagCompileTime::use_probe)
        {
          probe_at_t.set(properties.hydraulic_time);
        }

        KOKKOS_ASSERT(Kokkos::atomic_load(&current_container.n_cells) > 0);
        Kokkos::atomic_decrement(&current_container.n_cells);
        // clearState(MC::CellStatus::OUT);
        properties.reset();
        status = MC::CellStatus::OUT;
      }
    }
  }

  template <typename ListType, typename ResultViewType>
  KOKKOS_INLINE_FUNCTION size_t
  Kernel<ListType, ResultViewType>::__find_next_compartment(
      const std::size_t i_compartment,
      const double random_number,
      const Kokkos::View<const size_t *, Kokkos::LayoutStride, ComputeSpace>
          &i_neighbor) const
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor(0); // Default to the first neighbor

    // Iterate through the neighbors to find the appropriate next compartment
    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      // Get the cumulative probability range for the current neighbor
      const auto pi = cumulative_probability(i_compartment, k_neighbor);
      const auto pn = cumulative_probability(i_compartment, k_neighbor + 1);

      // Use of a Condition mask to avoid branching.
      // As the condition is not complex, ternary op manually guarantee that no
      // branching/warp divergence occurs
      next = (random_number <= pn && pi <= random_number)
                 ? i_neighbor(k_neighbor + 1)
                 : next;
    }

    return next; // Return the index of the chosen next compartment
  }

  //TODO CHECK IF WORTH IT 
  /*
  template <typename ListType, typename ResultViewType>
  KOKKOS_INLINE_FUNCTION size_t
  Kernel<ListType, ResultViewType>::_find_next_compartment_2(
      const std::size_t i_compartment,
      const double random_number,
      const NeighborsViewCompute
          &i_neighbor) const
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor(i_compartment,0); // Default to the first neighbor

    // Iterate through the neighbors to find the appropriate next compartment
    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      // Get the cumulative probability range for the current neighbor
      const auto pi = cumulative_probability(i_compartment, k_neighbor);
      const auto pn = cumulative_probability(i_compartment, k_neighbor + 1);

      // Use of a Condition mask to avoid branching.
      // As the condition is not complex, ternary op manually guarantee that no
      // branching/warp divergence occurs
      next = (random_number <= pn && pi <= random_number)
                 ? i_neighbor(i_compartment,k_neighbor + 1)
                 : next;
    }

    return next; // Return the index of the chosen next compartment
  }*/



} // namespace Simulation::KernelInline

#endif
