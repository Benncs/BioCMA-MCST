#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include "cma_read/light_2d_view.hpp"
#include "common/common.hpp"
#include "common/kokkos_vector.hpp"
#include <mc/particles/particle_model.hpp>
#include <cstdint>
#include <simulation/pc_hydro.hpp>
#include <span>
#include <simulation/alias.hpp>
// TODO REMOVE
#include <Kokkos_Core.hpp>
#include <cstddef>
#include <vector> 

namespace Simulation
{
  struct EigenData
  {
    Eigen::MatrixXd concentration;
    Eigen::MatrixXd total_mass;

    Eigen::MatrixXd biomass_contribution;

    Eigen::MatrixXd feed;
  };

  class ScalarSimulation
  {
  public:
    ScalarSimulation(ScalarSimulation &&other) noexcept;
    ScalarSimulation(const ScalarSimulation &other) noexcept = delete;
    ScalarSimulation(size_t n_compartments, size_t n_species, std::span<double> volume);
    ScalarSimulation operator=(const ScalarSimulation &other) = delete;
    ScalarSimulation operator=(ScalarSimulation &&other) = delete;

    ~ScalarSimulation() = default;

    bool deep_copy_concentration(const std::vector<double> &data);

    // Getters

    [[nodiscard]] auto get_device_concentration() const
    {
      return compute_concentration;
    }

    [[nodiscard]] std::span<double const> getVolumeData() const;

    [[nodiscard]] std::span<double> getContributionData();

    Eigen::DiagonalMatrix<double, -1> &getVolume();

    std::span<double> getConcentrationData();

    [[nodiscard]] Eigen::ArrayXXd getConcentrationArray() const;

    CmaRead::L2DView<double> getConcentrationView();

    [[nodiscard]] std::size_t n_species() const;

    [[nodiscard]] std::size_t n_row() const;
    [[nodiscard]] std::size_t n_col() const;

    [[nodiscard]] kernelContribution get_kernel_contribution() const;

    // Setters

    void set_mass();
    void set_kernel_contribs_to_host(kernelContribution c) const;

    void set_feed(std::uint64_t i_r, std::uint64_t i_c, double val);
    void set_sink(std::uint64_t i_compartment, double val);

    void setVolumes(std::span<const double> volumes, std::span<const double> inv_volumes);

    Eigen::ArrayXXd vec_kla; // TODO : Clean this

    Eigen::MatrixXd biomass_contribution;

    Eigen::MatrixXd mass_transfer;

    void performStep(double d_t, const FlowMatrixType &m_transition, const Eigen::MatrixXd &transfer_gas_liquid);

   

  private:
    Eigen::MatrixXd total_mass;

    Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace> host_concentration;
    Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace> compute_concentration;

    Eigen::DiagonalMatrix<double, -1> volumes_inverse;
    Eigen::DiagonalMatrix<double, -1> m_volumes;

    Eigen::SparseMatrix<double> feed;

    // Eigen::SparseMatrix<double> sink;

    Eigen::DiagonalMatrix<double, -1> sink;

    Eigen::MatrixXd alloc_concentrations;
    size_t n_r;
    size_t n_c;
    CmaRead::L2DView<double> view;

    Kokkos::View<double **, Kokkos::LayoutLeft, HostSpace> host_view_biomass_contribution;
  };

  inline size_t ScalarSimulation::n_col() const
  {
    return n_c;
  }

  inline Eigen::ArrayXXd ScalarSimulation::getConcentrationArray() const
  {
    return alloc_concentrations.array();
  }
  inline size_t ScalarSimulation::n_row() const
  {
    return n_r;
  }

  inline kernelContribution ScalarSimulation::get_kernel_contribution() const
  {
    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_view_biomass_contribution);
  }

  inline void ScalarSimulation::set_kernel_contribs_to_host(kernelContribution c) const
  {
    Kokkos::deep_copy(host_view_biomass_contribution, c);
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

  inline Eigen::DiagonalMatrix<double, -1> &ScalarSimulation::getVolume()
  {
    return m_volumes;
  }

  inline size_t ScalarSimulation::n_species() const
  {
    return this->n_r;
  }

  inline std::span<double> ScalarSimulation::getConcentrationData()
  {

    return {this->alloc_concentrations.data(), static_cast<size_t>(this->alloc_concentrations.size())};
  }

  inline std::span<double> ScalarSimulation::getContributionData()
  {
    return {this->biomass_contribution.data(), static_cast<size_t>(this->biomass_contribution.size())};
  }

  inline std::span<double const> ScalarSimulation::getVolumeData() const
  {
    return {m_volumes.diagonal().data(), static_cast<size_t>(m_volumes.rows())};
  }

  inline void ScalarSimulation::setVolumes(std::span<const double> volumes, std::span<const double> inv_volumes)
  {

    // SIGFAULT ?
    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(volumes.data(), static_cast<int>(volumes.size()));

    this->volumes_inverse.diagonal() = Eigen::Map<const Eigen::VectorXd>(inv_volumes.data(), static_cast<int>(inv_volumes.size()));
  }

  inline ScalarSimulation *makeScalarSimulation(size_t n_compartments, size_t n_species, std::span<double> volumes)
  {
    return new ScalarSimulation(n_compartments, n_species, volumes);
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__