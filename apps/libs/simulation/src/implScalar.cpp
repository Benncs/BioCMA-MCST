#include <Eigen/Dense>
#include <scalar_simulation.hpp>

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(ScalarSimulation &&other) noexcept
      : C(std::move(other.C)), Mtot(std::move(other.Mtot)), V(other.V)
  {
  }

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     size_t n_thread)
      : n_r(n_species), n_c(n_compartments)
  {
    Eigen::setNbThreads(static_cast<int>(n_thread));
    V = Eigen::DiagonalMatrix<double, -1>(static_cast<int>(n_compartments));
    V.setIdentity();

    this->Mtot = Eigen::MatrixXd(n_species, n_compartments);
    this->Mtot.setZero();

    this->C = Eigen::MatrixXd(n_species, n_compartments);
    this->C.setZero();

    this->vec_kla = Eigen::ArrayXXd(n_species,n_compartments);
    vec_kla.setZero();

    biomass_contribution = Eigen::MatrixXd(n_species, n_compartments);
    biomass_contribution.setZero();

    contribs.resize(n_thread);
    std::generate_n(contribs.begin(),
                    n_thread,
                    [this]()
                    {
                      auto m = Eigen::MatrixXd(C.rows(), C.cols());
                      m.setZero();
                      return m;
                    });
  }

  void ScalarSimulation::performStep(double d_t,
                                     const FlowMatrixType &m_transition,
                                     const Eigen::MatrixXd &transfer_gas_liquid)
  {

  
    Eigen::MatrixXd v_inverse = V.inverse();

    Mtot += d_t * (Mtot * v_inverse * m_transition + biomass_contribution + transfer_gas_liquid);

    C = Mtot * v_inverse;
  }
} // namespace Simulation
