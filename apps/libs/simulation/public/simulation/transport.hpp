#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__


#include <functional>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <simulation/matflows.hpp>


namespace Simulation
{

  FlowMatrixType FlowmapToMat(std::span<double> data, size_t n_row);

  FlowMatrixType get_transition_matrix(const FlowMatrixType &flows);

  Eigen::MatrixXd get_CP(const std::vector<std::vector<size_t>> &neighbors,
                         int nb_zone,
                         const Simulation::MatFlow &flows);

  using move_kernel =
      std::function<void(double, double, MC::Particles &, double)>;

  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer &container,
                  const MatFlow *flows);
} // namespace Simulation

#endif //__TRANSPORT_HPP__