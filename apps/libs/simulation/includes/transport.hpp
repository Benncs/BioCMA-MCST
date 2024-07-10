#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include "mc/domain.hpp"
#include <cma_read/flowmap.hpp>
#include <functional>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <simulation/pc_hydro.hpp>

namespace Simulation
{


  //  using move_kernel_t = std::function<void(
  //     double, double, MC::ReactorDomain &, MC::Particles &, double)>;

  // move_kernel_t population_balance_flow(MC::ReactorDomain &domain,
  //                                       const PreCalculatedHydroState *flows);

  FlowMatrixType
  get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows);

  Eigen::MatrixXd get_CP(CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                         const FlowMatrixType &m_transition);

  void kernel_exit(double d_t,
                   double random_number,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle);

  void kernel_move(double random_number,
                   double random_number2,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle,
                   double d_t,
                   const FlowMatrixType &m_transition,
                   const Eigen::MatrixXd &cumulative_probability);

} // namespace Simulation

#endif //__TRANSPORT_HPP__