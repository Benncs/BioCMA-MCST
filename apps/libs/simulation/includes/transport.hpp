#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include "cma_read/light_2d_view.hpp"
#include "common/thread_safe_container.hpp"
#include "mc/domain.hpp"
#include <Kokkos_Core.hpp>
#include <cma_read/flowmap.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <simulation/pc_hydro.hpp>

namespace
{
  inline bool probability_leaving(double random_number,
                                  double volume,
                                  double flow,
                                  double dt)
  {
    return (dt * flow / volume) > (-std::log(1 - random_number));
  }

  

} // namespace

namespace Simulation
{
  /**
   * @brief Finds the next compartment for a particle based on a random number
   * and cumulative probabilities.
   *
   * This function determines the next compartment for a particle by comparing a
   * random number against the cumulative probabilities of the neighboring
   * compartments. It returns the index of the chosen next compartment based on
   * the transition probabilities.
   *
   * @param i_compartment Index of the current compartment.
   * @param random_number Random number used to decide the next compartment.
   * @param i_neighbor Span of indices representing neighboring compartments.
   * @param cumulative_probability 2D view of cumulative probabilities for the
   * compartments.
   *
   * @return The index of the next compartment for the particle.
   */
  size_t
  find_next_compartment(int i_compartment,
                        double random_number,
                        std::span<const size_t> i_neighbor,
                        CmaRead::L2DView<const double> cumulative_probability);

  FlowMatrixType
  get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows);

  std::vector<double> get_diag_transition(const FlowMatrixType &m_transition);

  

  

  

} // namespace Simulation

#endif //__TRANSPORT_HPP__