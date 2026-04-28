#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include "mc/alias.hpp"
#include <common/eigen_diag.hpp>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
EIGEN_DIAG_POP

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <cstdint>
#include <eigen_kokkos.hpp>
#include <mc/traits.hpp>
#include <simulation/mass_transfer.hpp>
#include <span>
#include <vector>

using FlowMatrixType = Eigen::SparseMatrix<double>;

namespace Simulation
{

  class ScalarSimulation
  {
  public:
    ScalarSimulation(size_t n_compartments,
                     size_t n_species,
                     std::span<double> volume);

    ScalarSimulation(ScalarSimulation&& other) noexcept = delete;
    ScalarSimulation(const ScalarSimulation& other) noexcept = delete;
    ScalarSimulation operator=(const ScalarSimulation& other) = delete;
    ScalarSimulation operator=(ScalarSimulation&& other) = delete;
    ~ScalarSimulation() = default;

    bool deep_copy_concentration(const std::vector<double>& data);

    void reduce_contribs(std::span<const double> data);

    void set_transition(CmaUtils::StateCooMatrixType&& transition);

    void performStepGL(double d_t,
                       const ColMajorMatrixtype<double>& mtr,
                       MassTransfer::Sign sign);

    void performStep(double d_t);

    void synchro_sources();

    // void performStep(double d_t, const FlowMatrixType& m_transition);

    // Getters
    [[nodiscard]] ColMajorMatrixtype<double>& get_concentration() noexcept;
    [[nodiscard]] const DiagonalType& getVolume() const noexcept;
    [[nodiscard]] std::size_t n_row() const noexcept;
    [[nodiscard]] std::size_t n_col() const noexcept;

    [[nodiscard]] ColMajorKokkosScalarMatrix<double>
    get_device_concentration() const;

    [[nodiscard]] std::span<double const> volume_span() const;

    [[nodiscard]] std::span<double> contribution_span() const;

    [[nodiscard]] std::span<double> contribution_span_mut();

    [[nodiscard]] auto getConcentrationArray() const;
    [[nodiscard]] kernelContribution get_kernel_contribution() const;

    [[nodiscard]] std::span<double> getConcentrationData();

    [[nodiscard]] std::span<const double> getConcentrationData() const;

    // Setters

    void set_mass();
    void set_kernel_contribs_to_host() const;
    void set_feed(std::uint64_t i_r, std::uint64_t i_c, double val);
    void set_sink(std::uint64_t i_compartment, double val);
    void set_zero_contribs();
    void setVolumes(std::span<const double> volumes,
                    std::span<const double> inv_volumes);

    void clearNegs();

  private:
    std::size_t n_r;
    std::size_t n_c;
    ColMajorMatrixtype<double> total_mass;
    kernelContribution contribs;

    FlowMatrixType m_transition;

    DiagonalType volumes_inverse;
    DiagonalType m_volumes;
    DiagonalType sink;
    EigenKokkos<double> concentrations;
    RowMajorEigenKokkos<double> sources;

    // Kokkos::View<float**, ComputeSpace::array_layout, ComputeSpace> kc;
  };

  inline auto
  ScalarSimulation::getConcentrationArray() const
  {
    // return alloc_concentrations.array();
    return concentrations.eigen_data.array();
  }

  inline kernelContribution
  ScalarSimulation::get_kernel_contribution() const
  {
    // return sources.compute;
    return contribs;
  }

  inline void
  ScalarSimulation::set_zero_contribs()
  {
    sources.eigen_data.setZero();
    Kokkos::deep_copy(contribs, 0);
    Kokkos::deep_copy(sources.compute, 0);
    this->sink.setZero();
  }

  inline void
  ScalarSimulation::set_kernel_contribs_to_host() const
  {
    sources.update_compute_to_host();
  }

  inline void
  ScalarSimulation::set_feed(uint64_t i_r, uint64_t i_c, double val)
  {
    this->sources.eigen_data.coeffRef(EIGEN_INDEX(i_r), EIGEN_INDEX(i_c))
        += val;
  }

  inline void
  ScalarSimulation::set_sink(uint64_t i_compartment, double val)
  {
    this->sink.diagonal().coeffRef(EIGEN_INDEX(i_compartment)) += val;
  }

  inline const DiagonalType&
  ScalarSimulation::getVolume() const noexcept
  {
    return m_volumes;
  }

  inline std::span<double>
  ScalarSimulation::getConcentrationData()
  {
    return this->concentrations.get_span();
  }

  inline std::span<const double>
  ScalarSimulation::getConcentrationData() const
  {
    return this->concentrations.get_span();
  }

  inline std::span<double>
  ScalarSimulation::contribution_span() const
  {
    return { this->sources.host.data(),
             static_cast<size_t>(this->sources.host.size()) };
  }

  inline std::span<double>
  ScalarSimulation::contribution_span_mut()
  {
    return { this->sources.host.data(),
             static_cast<size_t>(this->sources.host.size()) };
  }

  inline std::span<double const>
  ScalarSimulation::volume_span() const
  {
    return { m_volumes.diagonal().data(),
             static_cast<size_t>(m_volumes.rows()) };
  }

  inline void
  ScalarSimulation::setVolumes(std::span<const double> volumes,
                               std::span<const double> inv_volumes)
  {
    KOKKOS_ASSERT(volumes.size() == inv_volumes.size()
                  && volumes.size() == n_col() && "scalar:setvolume")
    // SIGFAULT ?
    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    this->volumes_inverse.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        inv_volumes.data(), static_cast<int>(inv_volumes.size()));
  }

  inline ScalarSimulation*
  makeScalarSimulation(size_t n_compartments,
                       size_t n_species,
                       std::span<double> volumes)
  {
    return new ScalarSimulation(n_compartments, n_species, volumes); // NOLINT
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__
