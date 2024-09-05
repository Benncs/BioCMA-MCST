#ifndef __SIMULATION_MC_KERNEL_HPP
#define __SIMULATION_MC_KERNEL_HPP

#include <Kokkos_Core.hpp>
#include <cassert>
#include <common/kokkos_vector.hpp>
#include <mc/container_state.hpp>
#include <mc/events.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>

#ifndef NDEBUG
#  define ENABLE_KOKKOS_PROFILING
#  include <Kokkos_Profiling_ScopedRegion.hpp>
#endif

namespace Simulation::KernelInline
{
  KOKKOS_INLINE_FUNCTION bool probability_leaving(double random_number,
                                                  double volume,
                                                  double flow,
                                                  double dt)
  {
    // Use of Kokkos log to be sure that log is well defined
    return (dt * flow / volume) > (-Kokkos::log(1 - random_number));
  }

  KOKKOS_INLINE_FUNCTION size_t __find_next_compartment(
      const std::size_t i_compartment,
      double random_number,
      const Kokkos::View<const size_t *, Kokkos::LayoutStride, ComputeSpace>
          &i_neighbor,
      const Kokkos::View<const double **, Kokkos::LayoutStride, ComputeSpace>
          &cumulative_probability)
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor(0); // Default to the first neighbor

    // Iterate through the neighbors to find the appropriate next compartment
    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      // Get the cumulative probability range for the current neighbor
      const auto pi = cumulative_probability(i_compartment, k_neighbor);
      const auto pn = cumulative_probability(i_compartment, k_neighbor + 1);

      // Check if the random number falls within the probability range
      if (random_number <= pn && pi <= random_number)
      {
        next = i_neighbor(k_neighbor + 1); // Update to the next neighbor
        // No need to break, as we're looking for the last valid neighbor in the
        // range
      }
    }

    return next; // Return the index of the chosen next compartment
  }

  // TODO : CHECK THIS
  template <ParticleModel Model>
  KOKKOS_INLINE_FUNCTION void handle_exit(double random_number,
                                          double d_t,
                                          MC::Particle<Model> &particle,
                                          auto &&local_compartments,
                                          auto &&index,
                                          auto &&flow)
  {

    if (particle.properties.current_container != index ||
        particle.properties.status != MC::CellStatus::IDLE)
    {
      return;
    }

    auto &current_container = local_compartments(index);
    if (KernelInline::probability_leaving(
            random_number, current_container.volume_liq, flow, d_t))
    {
      particle.properties.status = MC::CellStatus::OUT;
      assert(Kokkos::atomic_load(&current_container.n_cells) > 0);
      Kokkos::atomic_decrement(&current_container.n_cells);
      particle.clearState(MC::CellStatus::OUT);
    }
  }


    template <ParticleModel Model>
  KOKKOS_INLINE_FUNCTION void handle_move(const std::size_t i_compartment,
                                          MC::Particle<Model> &particle,
                                          auto &local_compartments,
                                          auto &neighbors,
                                          auto &diag_transition,
                                          auto &cumulative_probability,
                                          auto &events,
                                          double d_t,auto&& local_rng)
  {

    auto &current_container = local_compartments(i_compartment);

    const auto i_neighbor =
        Kokkos::subview(neighbors, i_compartment, Kokkos::ALL);

    // Need 2 random numbers, use index 0 and 1 to acess to
    const auto random = local_rng.template double_uniform<2>();

    if (KernelInline::probability_leaving(random[0],
                                          current_container.volume_liq,
                                          diag_transition(i_compartment),
                                          d_t))
    {
      // Find the next compartment based on the random number and
      // cumulative
      // probabilities
      const size_t next = KernelInline::__find_next_compartment(
          i_compartment, random[1], i_neighbor, cumulative_probability);

      assert(Kokkos::atomic_load(&current_container.n_cells) > 0); // ??
      Kokkos::atomic_decrement(&current_container.n_cells);
      Kokkos::atomic_increment(&local_compartments(next).n_cells);
      particle.properties.current_container = next;
      events.template incr<MC::EventType::Move>();
    }
  };

} // namespace Simulation::KernelInline

#endif