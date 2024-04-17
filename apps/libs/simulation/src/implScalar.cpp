#include <common/common.hpp>
#include <Eigen/Dense>
#include <scalar_simulation.hpp>

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(ScalarSimulation &&other) noexcept
      : C(std::move(other.C)), Mtot(std::move(other.Mtot)), m_volumes(other.m_volumes),n_r(other.n_r),n_c(other.n_c)
  {
  }

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     size_t n_thread)
      : n_r(n_species), n_c(n_compartments)
  {
    Eigen::setNbThreads(EIGEN_INDEX(n_thread));
    m_volumes = Eigen::DiagonalMatrix<double, -1>(EIGEN_INDEX(n_compartments));
    m_volumes.setIdentity();

    volumes_inverse = Eigen::DiagonalMatrix<double, -1>(EIGEN_INDEX(n_compartments));
    volumes_inverse.setIdentity();

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

  
    Mtot += d_t * (Mtot * volumes_inverse * m_transition + biomass_contribution + transfer_gas_liquid);

    C = Mtot * volumes_inverse;
  }
} // namespace Simulation
