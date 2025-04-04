#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <common/common.hpp>

#include <scalar_simulation.hpp>
#include <simulation/alias.hpp>
#include <stdexcept>

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     std::span<double> volumes)
      : n_r(n_species), n_c(n_compartments), concentrations(n_r, n_c), sources(n_r, n_c)
  {
    if (volumes.size() != n_compartments)
    {
      throw std::invalid_argument("Volumes size mismatch");
    }

    const int n_row = EIGEN_INDEX(n_r);
    const int n_col = EIGEN_INDEX(n_c);

    m_volumes = Eigen::DiagonalMatrix<double, -1>(n_col);
    this->m_volumes.diagonal() =
        Eigen::Map<const Eigen::VectorXd>(volumes.data(), static_cast<int>(volumes.size()));

    volumes_inverse = Eigen::DiagonalMatrix<double, -1>(n_col);
    volumes_inverse.setIdentity();

    this->total_mass = ColMajorMatrixtype(n_row, n_col);
    this->total_mass.setZero();

    this->sink = DiagonalType(n_col);
    this->sink.setZero();
  }

  [[nodiscard]] ColMajorMatrixtype& ScalarSimulation::get_concentration()
  {
    return concentrations.eigen_data;
  }

  [[nodiscard]] ColMajorKokkosScalarMatrix ScalarSimulation::get_device_concentration() const
  {
    return concentrations.compute;
  }

  std::size_t ScalarSimulation::n_col() const
  {
    return n_c;
  }

  std::size_t ScalarSimulation::n_row() const
  {
    return n_r;
  }

  void ScalarSimulation::reduce_contribs(std::span<const double> data)
  {
    assert(data.size() == (n_c * n_r));
    sources.eigen_data.noalias() += Eigen::Map<Eigen::MatrixXd>(
        const_cast<double*>(data.data()), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
  }

  void ScalarSimulation::performStepGL(double d_t,
                                       const FlowMatrixType& m_transition,
                                       const ColMajorMatrixtype& mtr,
                                       MassTransfer::Sign sign)
  {
    PROFILE_SECTION("performStep_gl")
#define c concentrations.eigen_data

    total_mass = total_mass + d_t * (c * m_transition - c * sink + sources.eigen_data +
                                     static_cast<float>(sign) * mtr);
    c = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  void ScalarSimulation::performStep(double d_t, const FlowMatrixType& m_transition)
  {
    PROFILE_SECTION("performStep_gl")
#define c concentrations.eigen_data

    total_mass = total_mass + d_t * (c * m_transition - c * sink + sources.eigen_data);
    c = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  bool ScalarSimulation::deep_copy_concentration(const std::vector<double>& data)
  {
    if (data.size() != n_c * n_r)
    {
      std::cout<<data.size()<<std::endl;
      std::cout<<n_c * n_r<<std::endl;
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
