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
#include <simulation/alias.hpp>
#include <simulation/pc_hydro.hpp>
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

    [[nodiscard]] KokkosScalarMatrix<ComputeSpace> get_device_concentration() const;

    void reduce_contribs(std::span<const double> data);

    [[nodiscard]] std::span<double const> getVolumeData() const;

    [[nodiscard]] std::span<double> getContributionData() const;

    DiagonalType& getVolume();

    std::span<double> getConcentrationData();

    [[nodiscard]] Eigen::ArrayXXd getConcentrationArray() const;

    CmaRead::L2DView<double> getConcentrationView();

    [[nodiscard]] std::size_t n_species() const;

    [[nodiscard]] std::size_t n_row() const;
    [[nodiscard]] std::size_t n_col() const;

    [[nodiscard]] kernelContribution get_kernel_contribution() const;

    // Setters

    void set_mass();
    void set_kernel_contribs_to_host() const;

    void set_feed(std::uint64_t i_r, std::uint64_t i_c, double val);
    void set_sink(std::uint64_t i_compartment, double val);

    void setVolumes(std::span<const double> volumes, std::span<const double> inv_volumes);

    Eigen::ArrayXXd vec_kla; // TODO : Clean this

    void performStep(double d_t,
                     const FlowMatrixType& m_transition,
                     const MatrixType& transfer_gas_liquid);

    MatrixType& set_mass_transfer(MatrixType&& mtr);

    void set_zero_contribs()
    {
      b_contribs.eigen_data.setZero();
      Kokkos::deep_copy(b_contribs.compute, 0);
    }

    [[nodiscard]] const MatrixType& get_mass_transfer() const
    {
      return mass_transfer;
    };

  private:
    std::size_t n_r;
    std::size_t n_c;
    MatrixType total_mass;
    MatrixType mass_transfer;
    DiagonalType volumes_inverse;
    DiagonalType m_volumes;
    DiagonalType sink;
    SparseMatrixType feed;

    CmaRead::L2DView<double> view;
    EigenKokkos concentrations;
    EigenKokkos b_contribs;
  };

  inline size_t ScalarSimulation::n_col() const
  {
    return n_c;
  }

  inline Eigen::ArrayXXd ScalarSimulation::getConcentrationArray() const
  {
    // return alloc_concentrations.array();
    return concentrations.eigen_data.array();
  }
  inline size_t ScalarSimulation::n_row() const
  {
    return n_r;
  }

  inline kernelContribution ScalarSimulation::get_kernel_contribution() const
  {
    return b_contribs.compute;
  }

  inline void ScalarSimulation::set_kernel_contribs_to_host() const
  {

    b_contribs.update_compute_to_host();

    // Kokkos::deep_copy(host_view_biomass_contribution, kernel_view_biomass_contribution);
    // Kokkos::deep_copy(kernel_view_biomass_contribution, 0);
  }

  inline void ScalarSimulation::set_feed(uint64_t i_r, uint64_t i_c, double val)
  {
    this->feed.coeffRef(EIGEN_INDEX(i_r), EIGEN_INDEX(i_c)) = val;
  }

  inline void ScalarSimulation::set_sink(uint64_t i_compartment, double val)
  {
    this->sink.diagonal().coeffRef(EIGEN_INDEX(i_compartment)) = val;
  }

  inline CmaRead::L2DView<double> ScalarSimulation::getConcentrationView()
  {
    return view;
  }

  inline Eigen::DiagonalMatrix<double, -1>& ScalarSimulation::getVolume()
  {
    return m_volumes;
  }

  inline size_t ScalarSimulation::n_species() const
  {
    return this->n_r;
  }

  inline std::span<double> ScalarSimulation::getConcentrationData()
  {
    return this->concentrations.get_span();
    // return {this->alloc_concentrations.data(),
    //         static_cast<size_t>(this->alloc_concentrations.size())};
  }

  inline std::span<double> ScalarSimulation::getContributionData() const
  {
    // return {this->biomass_contribution.data(),
    //         static_cast<size_t>(this->biomass_contribution.size())};
    return {this->b_contribs.host.data(), static_cast<size_t>(this->b_contribs.host.size())};
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