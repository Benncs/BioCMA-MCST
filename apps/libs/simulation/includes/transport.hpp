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

  using move_kernel = std::function<void(
      double, double, MC::ReactorDomain &, MC::Particles &, double)>;

  move_kernel population_balance_flow(MC::ReactorDomain &domain,
                                      MC::ParticlesContainer &container,
                                      const MatFlow *flows);
} // namespace Simulation

#endif //__TRANSPORT_HPP__