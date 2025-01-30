#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <scalar_simulation.hpp>
#include <simulation/alias.hpp>

EigenKokkos::EigenKokkos(std::size_t n_row, std::size_t n_col)
{
  eigen_data = Eigen::MatrixXd(n_row, n_col);
  eigen_data.setZero();

  host = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>(eigen_data.data(), n_row, n_col);
  compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host);
}

void EigenKokkos::update_host_to_compute() const
{
  Kokkos::deep_copy(compute, host);
}

void EigenKokkos::update_compute_to_host() const
{
  Kokkos::deep_copy(host, compute);
}

namespace Simulation
{

  // ScalarSimulation::ScalarSimulation(ScalarSimulation&& other) noexcept
  //     : total_mass(std::move(other.total_mass)), m_volumes(other.m_volumes),
  //       alloc_concentrations(std::move(other.alloc_concentrations)), n_r(other.n_r),
  //       n_c(other.n_c)
  // {
  // }

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     std::span<double> volumes)
      : n_r(n_species), n_c(n_compartments), concentrations(n_r, n_c), b_contribs(n_r, n_c)
  {

    const int n_row = EIGEN_INDEX(n_r);
    const int n_col = EIGEN_INDEX(n_c);

    m_volumes = Eigen::DiagonalMatrix<double, -1>(n_col);
    this->m_volumes.diagonal() =
        Eigen::Map<const Eigen::VectorXd>(volumes.data(), static_cast<int>(volumes.size()));

    volumes_inverse = Eigen::DiagonalMatrix<double, -1>(n_col);
    volumes_inverse.setIdentity();

    this->total_mass = Eigen::MatrixXd(n_row, n_col);
    this->total_mass.setZero();

    this->feed = Eigen::SparseMatrix<double>(n_row, n_col);
    this->feed.setZero();

    this->sink = Eigen::DiagonalMatrix<double, -1>(n_col);
    this->sink.setZero();

    this->vec_kla = Eigen::ArrayXXd(n_row, n_col);
    this->vec_kla.setZero();

    this->mass_transfer = Eigen::MatrixXd(n_row, n_col);
    this->mass_transfer.setZero();

    // todo remove
    view = CmaRead::L2DView<double>({this->concentrations.eigen_data.data(),
                                     static_cast<size_t>(this->concentrations.eigen_data.size())},
                                    this->concentrations.eigen_data.rows(),
                                    this->concentrations.eigen_data.cols(),
                                    false);
  }

  [[nodiscard]] KokkosScalarMatrix<ComputeSpace> ScalarSimulation::get_device_concentration() const
  {
    return concentrations.compute;
    // return compute_concentration;
  }

  void ScalarSimulation::reduce_contribs(std::span<const double> data)
  {
    assert(data.size() == (n_c * n_r));
    b_contribs.eigen_data.noalias() += Eigen::Map<Eigen::MatrixXd>(
        const_cast<double*>(data.data()), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
  }

  void ScalarSimulation::performStep(double d_t,
                                     const FlowMatrixType& m_transition,
                                     const Eigen::MatrixXd& transfer_gas_liquid)
  {
    PROFILE_SECTION("performStep")
    auto& c = concentrations.eigen_data;
    total_mass +=
        d_t * (c * m_transition - c * sink + b_contribs.eigen_data + feed + transfer_gas_liquid);

    // total_mass = (total_mass.array()<0.).select(0,total_mass);
    // total_mass += d_t * (alloc_concentrations * m_transition -
    //                      alloc_concentrations * sink + feed);

    c = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  Eigen::MatrixXd& ScalarSimulation::set_mass_transfer(Eigen::MatrixXd&& mtr)
  {
    this->mass_transfer = std::move(mtr);
    return this->mass_transfer;
  }

  bool ScalarSimulation::deep_copy_concentration(const std::vector<double>& data)
  {
    if (data.size() != n_c * n_r)
    {
      return false;
    }

    Eigen::Map<const Eigen::MatrixXd> temp_map(data.data(), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
    this->concentrations.eigen_data = temp_map; // Performs deep copy
    return true;
  }

  void ScalarSimulation::set_mass()
  {
    total_mass = this->concentrations.eigen_data * m_volumes;
  }
} // namespace Simulation
