#include <Eigen/Dense>
#include <simulation/scalar_simulation.hpp>

namespace Simulation
{
  ScalarSimulation::ScalarSimulation(ScalarSimulation &&other) noexcept
  {
    this->Mtot = std::move(other.Mtot);
    this->V = std::move(other.V);
    this->C = std::move(other.C);
  }

  ScalarSimulation::ScalarSimulation(size_t n_compartments, size_t n_species)
  { 

    V = Eigen::DiagonalMatrix<double,-1>(n_compartments);   //FIXEME this should be diag or sparse 
    V.setIdentity();

    this->Mtot = Eigen::MatrixXd(n_species, n_compartments);
    this->Mtot.setZero();

    this->C = Eigen::MatrixXd(n_species, n_compartments);
    this->C.setZero();

    for (size_t i = 0; i < n_compartments; ++i)
    {
      C.coeffRef(0, static_cast<int>(i)) = static_cast<double>(i);
    }

    Mtot = C*V;
  }

  void
  ScalarSimulation::performStep(double d_t,
                                const Eigen::SparseMatrix<double> &m_transition,
                                const Eigen::MatrixXd &phi_v)
  {
    // self.C = self.C+ dt*np.dot(self.C,m_transition) +
    // np.dot(np.sum(self.phi_v,axis=0),self._mat_v)
    Mtot = Mtot + d_t * (C * m_transition) + phi_v*V;
  }
} // namespace Simulation