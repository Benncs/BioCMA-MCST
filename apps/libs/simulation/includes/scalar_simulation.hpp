#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <cma_read/light_2d_view.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/particles/particle_model.hpp>
#include <pc_hydro.hpp>
#include <simulation/alias.hpp>
#include <span>
#include <vector>

constexpr auto DataLayoutEigen = Eigen::ColMajor;
constexpr auto CompileMatrixSizeEigen = -1;
using MatrixType =
    Eigen::Matrix<double, CompileMatrixSizeEigen, CompileMatrixSizeEigen, DataLayoutEigen>;
using SparseMatrixType = Eigen::SparseMatrix<double, DataLayoutEigen>;
using DiagonalType = Eigen::DiagonalMatrix<double, CompileMatrixSizeEigen>;

template <typename ExecSpace>
using KokkosScalarMatrix = Kokkos::View<double**, Kokkos::LayoutLeft, ExecSpace>;

struct EigenKokkos
{
  KokkosScalarMatrix<HostSpace> host;
  KokkosScalarMatrix<ComputeSpace> compute;
  MatrixType eigen_data;

  EigenKokkos(std::size_t n_row, std::size_t n_col);
  [[nodiscard]] std::span<const double> get_span() const
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }
  std::span<double> get_span()
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }

  void update_host_to_compute() const;
  void update_compute_to_host() const;
};

namespace Simulation
{

  class ScalarSimulation
  {
  public:
    ScalarSimulation(size_t n_compartments, size_t n_species, std::span<double> volume);

    ScalarSimulation(ScalarSimulation&& other) noexcept = delete;
    ScalarSimulation(const ScalarSimulation& other) noexcept = delete;
    ScalarSimulation operator=(const ScalarSimulation& other) = delete;
    ScalarSimulation operator=(ScalarSimulation&& other) = delete;
    ~ScalarSimulation() = default;

    bool deep_copy_concentration(const std::vector<double>& data);
    void reduce_contribs(std::span<const double> data);

    void performStep(double d_t,
                     const FlowMatrixType& m_transition,
                     const MatrixType& transfer_gas_liquid);

    void performStep(double d_t, const FlowMatrixType& m_transition);

    // Getters
    [[nodiscard]] KokkosScalarMatrix<ComputeSpace> get_device_concentration() const;
    [[nodiscard]] std::span<double const> getVolumeData() const;
    [[nodiscard]] std::span<double> getContributionData() const;
    [[nodiscard]] const DiagonalType& getVolume() const;
    [[nodiscard]] Eigen::ArrayXXd getConcentrationArray() const;
    [[deprecated("remove depenency to cmaread")]] CmaRead::L2DView<double> getConcentrationView();
    [[nodiscard]] kernelContribution get_kernel_contribution() const;
    [[nodiscard]] const MatrixType& get_mass_transfer() const;
    [[nodiscard]] std::span<double> getConcentrationData();
    

    [[nodiscard]] std::size_t n_row() const;
    [[nodiscard]] std::size_t n_col() const;

    // Setters

    void set_mass();
    void set_kernel_contribs_to_host() const;
    void set_feed(std::uint64_t i_r, std::uint64_t i_c, double val);
    void set_sink(std::uint64_t i_compartment, double val);
    void set_zero_contribs();
    void setVolumes(std::span<const double> volumes, std::span<const double> inv_volumes);

  private:
    std::size_t n_r;
    std::size_t n_c;
    MatrixType total_mass;
    
    DiagonalType volumes_inverse;
    DiagonalType m_volumes;
    DiagonalType sink;

    CmaRead::L2DView<double> view;
    EigenKokkos concentrations;
    EigenKokkos sources;

    
  };

  inline Eigen::ArrayXXd ScalarSimulation::getConcentrationArray() const
  {
    // return alloc_concentrations.array();
    return concentrations.eigen_data.array();
  }

  inline kernelContribution ScalarSimulation::get_kernel_contribution() const
  {
    return sources.compute;
  }

  inline void ScalarSimulation::set_zero_contribs()
  {
    sources.eigen_data.setZero();
    Kokkos::deep_copy(sources.compute, 0);
  }



  inline void ScalarSimulation::set_kernel_contribs_to_host() const
  {
    sources.update_compute_to_host();
  }

  inline void ScalarSimulation::set_feed(uint64_t i_r, uint64_t i_c, double val)
  {
    this->sources.eigen_data.coeffRef(EIGEN_INDEX(i_r), EIGEN_INDEX(i_c)) += val;
  }

  inline void ScalarSimulation::set_sink(uint64_t i_compartment, double val)
  {
    this->sink.diagonal().coeffRef(EIGEN_INDEX(i_compartment)) = val;
  }

  inline const DiagonalType& ScalarSimulation::getVolume() const
  {
    return m_volumes;
  }

  inline std::span<double> ScalarSimulation::getConcentrationData()
  {
    return this->concentrations.get_span();
  }

  inline std::span<double> ScalarSimulation::getContributionData() const
  {
    return {this->sources.host.data(), static_cast<size_t>(this->sources.host.size())};
  }

  inline std::span<double const> ScalarSimulation::getVolumeData() const
  {
    return {m_volumes.diagonal().data(), static_cast<size_t>(m_volumes.rows())};
  }

  inline void ScalarSimulation::setVolumes(std::span<const double> volumes,
                                           std::span<const double> inv_volumes)
  {

    // SIGFAULT ?
    this->m_volumes.diagonal() =
        Eigen::Map<const Eigen::VectorXd>(volumes.data(), static_cast<int>(volumes.size()));

    this->volumes_inverse.diagonal() =
        Eigen::Map<const Eigen::VectorXd>(inv_volumes.data(), static_cast<int>(inv_volumes.size()));
  }

  inline ScalarSimulation*
  makeScalarSimulation(size_t n_compartments, size_t n_species, std::span<double> volumes)
  {
    return new ScalarSimulation(n_compartments, n_species, volumes); // NOLINT
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__