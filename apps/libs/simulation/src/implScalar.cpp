#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <scalar_simulation.hpp>

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(ScalarSimulation &&other) noexcept
      : alloc_concentrations(std::move(other.alloc_concentrations)),
        total_mass(std::move(other.total_mass)), m_volumes(other.m_volumes),
        n_r(other.n_r), n_c(other.n_c)
  {
  }

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     std::span<double> volumes)
      : n_r(n_species), n_c(n_compartments)
  {

    const int n_row = EIGEN_INDEX(n_r);
    const int n_col = EIGEN_INDEX(n_c);

    m_volumes = Eigen::DiagonalMatrix<double, -1>(n_col);
    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    volumes_inverse = Eigen::DiagonalMatrix<double, -1>(n_col);
    volumes_inverse.setIdentity();

    this->total_mass = Eigen::MatrixXd(n_row, n_col);
    this->total_mass.setZero();

    this->alloc_concentrations = Eigen::MatrixXd(n_row, n_col);
    this->alloc_concentrations.setZero();

    this->feed = Eigen::SparseMatrix<double>(n_row, n_col);
    this->feed.setZero();

    this->sink = Eigen::DiagonalMatrix<double, -1>(n_col);
    this->sink.setZero();

    this->vec_kla = Eigen::ArrayXXd(n_row, n_col);
    this->vec_kla.setZero();

    this->mass_transfer = Eigen::MatrixXd(n_row, n_col);
    this->mass_transfer.setZero();

    biomass_contribution = Eigen::MatrixXd(n_row, n_col);
    biomass_contribution.setZero();

    view = CmaRead::L2DView<double>(
        {this->alloc_concentrations.data(),
         static_cast<size_t>(this->alloc_concentrations.size())},
        alloc_concentrations.rows(),
        alloc_concentrations.cols(),
        false);

    host_view_biomass_contribution =
        Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace>(
            biomass_contribution.data(),
            biomass_contribution.rows(),
            biomass_contribution.cols());

    host_concentration = Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace>(
        alloc_concentrations.data(),
        alloc_concentrations.rows(),
        alloc_concentrations.cols());

    compute_concentration =
        Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_concentration);
  }

  void ScalarSimulation::performStep(double d_t,
                                     const FlowMatrixType &m_transition,
                                     const Eigen::MatrixXd &transfer_gas_liquid)
  {    

    total_mass +=
        d_t * ( alloc_concentrations * m_transition - alloc_concentrations *sink + biomass_contribution + feed +
               transfer_gas_liquid );

    // total_mass += d_t * (alloc_concentrations * m_transition -
    //                      alloc_concentrations * sink + feed);

    alloc_concentrations = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    Kokkos::deep_copy(compute_concentration, host_concentration);
  }

  bool
  ScalarSimulation::deep_copy_concentration(const std::vector<double> &data)
  {
    if (data.size() != n_c * n_r)
    {
      return false;
    }

    Eigen::Map<const Eigen::MatrixXd> temp_map(
        data.data(), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
    this->alloc_concentrations = temp_map; // Performs deep copy

    return true;
  }

  void ScalarSimulation::set_mass()
  {
    total_mass = alloc_concentrations * m_volumes;
  }
} // namespace Simulation
