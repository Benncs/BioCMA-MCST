#ifndef __GET_CUMULATIVE_PROBABILITY_HPP__
#define __GET_CUMULATIVE_PROBABILITY_HPP__

#include <simulation/pc_hydro.hpp>

namespace Simulation
{

  /**
   * @brief Computes the cumulative probability matrix from the transition
   * matrix and neighbor information.
   *
   * This function calculates the cumulative probability matrix for each zone
   * (or compartment) based on the transition probabilities and the neighboring
   * information. The cumulative probability matrix is used for determining the
   * likelihood of moving to different compartments.
   *
   * @param neighbors A view of the neighbors for each compartment.
   * @param m_transition The transition matrix that provides the probabilities
   * of transitioning from one compartment to another.
   *
   * @return An Eigen::MatrixXd representing the cumulative probability matrix.
   */
  Eigen::MatrixXd get_cumulative_probabilities(
      CmaRead::Neighbors::Neighbors_const_view_t neighbors,
      const FlowMatrixType &m_transition);

} // namespace Simulation
#endif //__GET_CUMULATIVE_PROBABILITY_HPP__