#include <Eigen/Dense>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
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

    this->feed = Eigen::SparseMatrix<double>(n_species, n_compartments);
    this->feed.setZero();

    this->sink = Eigen::SparseMatrix<double>(n_compartments, n_compartments);
    this->sink.setZero();

    this->vec_kla = Eigen::ArrayXXd(n_species, n_compartments);
    this->vec_kla.setZero();

    biomass_contribution = Eigen::MatrixXd(n_species, n_compartments);
    biomass_contribution.setZero();

    view = CmaRead::L2DView<double>(
        {this->concentration.data(),
         static_cast<size_t>(this->concentration.size())},
        concentration.rows(),
        concentration.cols(),
        false);

    k_contribs = Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace>(
        biomass_contribution.data(),
        biomass_contribution.rows(),
        biomass_contribution.cols());

    host_concentration = Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace>(
        concentration.data(), concentration.rows(), concentration.cols());

    compute_concentration =
        Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_concentration);
  }

  void ScalarSimulation::performStep(double d_t,
                                     const FlowMatrixType &m_transition,
                                     const Eigen::MatrixXd &transfer_gas_liquid)
  {

    total_mass.noalias() +=
        d_t * (concentration * m_transition + biomass_contribution + feed + transfer_gas_liquid*m_volumes - concentration * sink);

    concentration = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    Kokkos::deep_copy(compute_concentration, host_concentration);
  }

  bool ScalarSimulation::deep_copy_liquid_concentration(
      const std::vector<double> &data)
  {
    if (data.size() != n_c * n_r)
    {
      return false;
    }

    Eigen::Map<const Eigen::MatrixXd> temp_map(
        data.data(), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
    this->concentration = temp_map; // Performs deep copy

    return true;
  }
} // namespace Simulation
