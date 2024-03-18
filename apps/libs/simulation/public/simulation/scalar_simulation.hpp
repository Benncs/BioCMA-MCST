#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Simulation
{
  class ScalarSimulation
  {
  public:
    ScalarSimulation(ScalarSimulation &&other) noexcept;
    ScalarSimulation(size_t n_compartments, size_t n_species);
    Eigen::MatrixXd C;
    Eigen::MatrixXd Mtot;
    Eigen::DiagonalMatrix<double, -1> V;
    void performStep(double d_t,
                     const Eigen::SparseMatrix<double> &m_transition,
                     const Eigen::MatrixXd &phi_v);

  private:
    void updateC();
    size_t n_r;
    size_t n_c;
  };
} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__