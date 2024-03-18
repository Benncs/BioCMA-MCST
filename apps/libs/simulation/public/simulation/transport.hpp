#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include <Eigen/Sparse>
#include <functional>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>

namespace Simulation
{
  Eigen::SparseMatrix<double>
  get_transition_matrix(const Eigen::SparseMatrix<double> &flows);

  struct MatFlow
  {
    Eigen::SparseMatrix<double> flows;
    Eigen::SparseMatrix<double> transition_matrix;
  };

  using move_kernel =
      std::function<void(double, double, MC::Particles &, double)>;

  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer &container,
                  const MatFlow &flows);
} // namespace Simulation

#endif //__TRANSPORT_HPP__