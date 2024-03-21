#include <Eigen/Dense>
#include <simulation/scalar_simulation.hpp>

#include <iostream>

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
    
    V = Eigen::DiagonalMatrix<double, -1>(static_cast<int>(n_compartments));
    V.setIdentity();

    this->Mtot = Eigen::MatrixXd(n_species, n_compartments);
    this->Mtot.setZero();

    this->C = Eigen::MatrixXd(n_species, n_compartments);
    this->C.setZero();

    biomass_contribution = Eigen::MatrixXd(n_species, n_compartments);
    biomass_contribution.setZero();

    for (size_t i = 0; i < n_compartments; ++i)
    {
      C.coeffRef(0, static_cast<int>(i)) = 150;
    }

    Mtot = C * V;
  }

  void
  ScalarSimulation::performStep(double d_t,
                                const Eigen::SparseMatrix<double> &m_transition)
  {
    
    auto& Vinv = V.inverse();

    Mtot = Mtot + d_t * (Mtot * Vinv * m_transition) + biomass_contribution*V;


    C = Mtot * Vinv;
  }
} // namespace Simulation
