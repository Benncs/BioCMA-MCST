#ifndef __CMA_UTILS_PROXY_CACHE_HPP__
#define __CMA_UTILS_PROXY_CACHE_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/reactorstate.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <common/kokkos_vector.hpp>
#include <vector>

namespace CmaUtils
{
  /**
   * @brief Proxy class for PreCalculatedHydroState structure, which is used to fill this struct
   * with data from an external library.
   */
  class ProxyPreCalculatedHydroState
  {
  public:
    /**
     * @brief Set the transition matrix using a constant flow map view.
     * @param flows_view Constant view of FlowMap from CmaRead.
     */
    void set_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows_view);

    /**
     * @brief Set the transition matrix using an rvalue reference to FlowMatrixType.
     * @param matrix Rvalue reference to a FlowMatrixType.
     */
    void set_transition_matrix(FlowMatrixType&& matrix);

    /**
     * @brief Set the cumulative probability using a constant neighbors view.
     * @param neighbors Constant view of Neighbors from CmaRead.
     */
    void set_cumulative_probability(const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);

    /**
     * @brief Set the diagonal elements of the transition matrix.
     * @param diag Rvalue reference to a vector of double containing diagonal elements.
     */
    void set_diag_transition(std::vector<double>&& diag);

    PreCalculatedHydroState
        state; ///< Instance of PreCalculatedHydroState containing the filled data.
  };
} // namespace CmaUtils

#endif