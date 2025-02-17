#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include <Kokkos_Core.hpp>
#include <cma_read/flowmap.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <cstddef>
#include <span>
#include <vector>

namespace CmaUtils
{
  FlowMatrixType get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows);

  std::vector<double> get_diag_transition(const FlowMatrixType& m_transition);

} // namespace CmaUtils

#endif //__TRANSPORT_HPP__