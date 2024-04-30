#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include "mc/domain.hpp"
#include <functional>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <simulation/matflows.hpp>

namespace Simulation
{

  FlowMatrixType flowmap_to_matrix(std::span<double> data, size_t n_row);

  FlowMatrixType get_transition_matrix(const FlowMatrixType &flows);

  Eigen::MatrixXd get_CP(const std::vector<std::vector<size_t>> &neighbors,
                         int nb_zone,
                         const Simulation::MatFlow &flows);

  using move_kernel_t = std::function<void(
      double, double, MC::ReactorDomain &, MC::Particles &, double)>;

  move_kernel_t population_balance_flow(MC::ReactorDomain &domain,
                                        const MatFlow *flows);

  void kernel_move(double random_number,
                   double random_number2,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle,
                   double d_t,
                   const FlowMatrixType &m_transition,
                   const Eigen::MatrixXd &cumulative_probability);

} // namespace Simulation

#endif //__TRANSPORT_HPP__