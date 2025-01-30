#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include <cma_read/light_2d_view.hpp>

#include <Kokkos_Core.hpp>
#include <cma_read/flowmap.hpp>
#include <cstddef>
#include <pc_hydro.hpp>
#include <span>
#include <vector>
namespace Simulation
{

  FlowMatrixType get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows);

  std::vector<double> get_diag_transition(const FlowMatrixType& m_transition);

} // namespace Simulation

#endif //__TRANSPORT_HPP__