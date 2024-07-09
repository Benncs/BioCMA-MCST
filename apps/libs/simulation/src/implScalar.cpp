#include <Eigen/Dense>
#include <common/common.hpp>
#include <scalar_simulation.hpp>

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(ScalarSimulation &&other) noexcept
      : concentration(std::move(other.concentration)),
        total_mass(std::move(other.total_mass)), m_volumes(other.m_volumes),
        n_r(other.n_r), n_c(other.n_c)
  {
  }

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     size_t n_thread,
                                     std::span<double> volumes)
      : n_r(n_species), n_c(n_compartments)
  {

    m_volumes = Eigen::DiagonalMatrix<double, -1>(EIGEN_INDEX(n_compartments));
    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    volumes_inverse =
        Eigen::DiagonalMatrix<double, -1>(EIGEN_INDEX(n_compartments));
    volumes_inverse.setIdentity();

    this->total_mass = Eigen::MatrixXd(n_species, n_compartments);
    this->total_mass.setZero();

    this->concentration = Eigen::MatrixXd(n_species, n_compartments);
    this->concentration.setZero();

    this->vec_kla = Eigen::ArrayXXd(n_species, n_compartments);
    vec_kla.setZero();

    biomass_contribution = Eigen::MatrixXd(n_species, n_compartments);
    biomass_contribution.setZero();

    contribs.resize(n_thread);
    std::generate_n(contribs.begin(),
                    n_thread,
                    [this]()
                    {
                      auto m = Eigen::MatrixXd(concentration.rows(),
                                               concentration.cols());
                      m.setZero();
                      return m;
                    });

    view = CmaRead::L2DView<double>(
        {this->concentration.data(),
         static_cast<size_t>(this->concentration.size())},
        concentration.rows(),
        concentration.cols(),
        false);
  }

  void ScalarSimulation::performStep(double d_t,
                                     const FlowMatrixType &m_transition,
                                     const Eigen::MatrixXd &transfer_gas_liquid,const Eigen::MatrixXd &feed)
  {

    total_mass.noalias() +=
        d_t * (concentration * m_transition + biomass_contribution + feed+
               (transfer_gas_liquid)*m_volumes);



    concentration.noalias() = total_mass * volumes_inverse;
  }

    
} // namespace Simulation
