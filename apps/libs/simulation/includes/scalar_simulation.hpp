#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include "Kokkos_Assert.hpp"
#include <simulation/mass_transfer.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/particles/particle_model.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <simulation/alias.hpp>
#include <span>
#include <vector>
#include <eigen_kokkos.hpp>
#include <mc/traits.hpp>

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

    void performStepGL(double d_t,
                     const FlowMatrixType& m_transition,
                     const MatrixType& mtr,MassTransfer::Sign sign);

    void performStep(double d_t,
                     const FlowMatrixType& m_transition);

    // void performStep(double d_t, const FlowMatrixType& m_transition);

    // Getters
    [[nodiscard]] MatrixType& get_concentration();
    [[nodiscard]] KokkosScalarMatrix<ComputeSpace> get_device_concentration() const;
    [[nodiscard]] std::span<double const> getVolumeData() const;
    [[nodiscard]] std::span<double> getContributionData() const;
    [[nodiscard]] const DiagonalType& getVolume() const;
    [[nodiscard]] auto getConcentrationArray() const;
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
    EigenKokkos concentrations;
    EigenKokkos sources;

    
  };

  inline auto ScalarSimulation::getConcentrationArray() const
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
    KOKKOS_ASSERT(volumes.size()==inv_volumes.size() && volumes.size()==n_col()&&"scalar:setvolume")
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