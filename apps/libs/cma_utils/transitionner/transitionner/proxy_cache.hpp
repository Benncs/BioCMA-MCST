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
  class ProxyPreCalculatedHydroState
  {
  public:
    void set_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows_view);
    void set_transition_matrix(FlowMatrixType&& matrix);
    void set_cumulative_probability(const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);
    void set_diag_transition(std::vector<double>&& diag);

    PreCalculatedHydroState state;
  };
} // namespace CmaUtils

#endif