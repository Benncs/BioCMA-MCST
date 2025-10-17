#ifndef __GET_CUMULATIVE_PROBABILITY_HPP__
#define __GET_CUMULATIVE_PROBABILITY_HPP__

#ifndef NDEBUG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wnan-infinity-disabled"
#endif 
#include <Eigen/Core>
#ifndef NDEBUG
#pragma GCC diagnostic pop
#endif 
#include <Kokkos_Core.hpp>
#include <cma_read/flowmap.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/neighbors.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <vector>

namespace CmaUtils
{
  using CumulativeProbaType = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
  /**
   * @brief Computes the cumulative probability matrix from the transition
   * matrix and neighbor information.
   *
   * This function calculates the cumulative probability matrix for each zone
   * (or compartment) based on the transition probabilities and the neighboring
   * information. The cumulative probability matrix is used for determining the
   * likelihood of moving to different compartments. This is by construction non-strictly increasing
   * over rows
   *
   * @param neighbors A view of the neighbors for each compartment.
   * @param m_transition The transition matrix that provides the probabilities
   * of transitioning from one compartment to another.
   *
   * @return An Eigen::MatrixXd representing the cumulative probability matrix.
   */
  CumulativeProbaType
  get_cumulative_probabilities(CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                               const FlowMatrixType& m_transition);

  /**
   * @brief Calculate the transition matrix representing phase flow.
   *
   * Converts a flow map read from the CmaRead private library into FlowMatrixType, a public
   * typedef.
   * @param flows Constant view of FlowMap from CmaRead.
   * @return FlowMatrixType representing phase flow.
   */
  FlowMatrixType get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows);

  /**
   * @brief Return the diagonal of the transition matrix (outflows).
   *
   * This function performs a copy of the diagonal elements.
   * @param m_transition Reference to the transition matrix.
   * @return Vector of double containing the diagonal elements.
   */
  std::vector<double> get_diag_transition(const FlowMatrixType& m_transition);

} // namespace CmaUtils
#endif //__GET_CUMULATIVE_PROBABILITY_HPP__