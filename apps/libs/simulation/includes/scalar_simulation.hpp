#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include "common/traits.hpp"
#include "mc/alias.hpp"
#include <common/eigen_diag.hpp>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <kokkos_eigen.hpp>
EIGEN_DIAG_POP

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/traits.hpp>
#include <simulation/mass_transfer.hpp>
#include <span>
#include <vector>

namespace Simulation
{
  template <FloatingPointType ftype>
  using FlowMatrixType = Eigen::SparseMatrix<ftype>;
  using mass_balance_float_type = double;
  class ScalarSimulation
  {
    using concentration_t
        = KokkosEigen::KokkosEigen2D<double,
                                     Kokkos::LayoutLeft,
                                     Kokkos::DefaultExecutionSpace>;

  public:
    using concentration_view_t = concentration_t::host_view_type;
    using concentration_float_type = concentration_view_t::value_type;

    ScalarSimulation(size_t n_compartments,
                     size_t n_species,
                     std::span<mass_balance_float_type> volume);

    ScalarSimulation(ScalarSimulation&& other) noexcept = delete;
    ScalarSimulation(const ScalarSimulation& other) noexcept = delete;
    ScalarSimulation operator=(const ScalarSimulation& other) = delete;
    ScalarSimulation operator=(ScalarSimulation&& other) = delete;
    ~ScalarSimulation() = default;

    bool
    deep_copy_concentration(const std::vector<concentration_float_type>& data);

    // void reduce_contribs(std::span<const mass_balance_float_type> data);

    void set_transition(CmaUtils::StateCooMatrixType&& transition);

    void performStepGL(
        double d_t,
        const KokkosEigen::Alias::ColMajorMatrixtype<mass_balance_float_type>&
            mtr,
        MassTransfer::Sign sign);

    void performStep(double d_t);

    void synchro_sources();

    // Getters

    [[nodiscard]] std::size_t n_row() const noexcept;
    [[nodiscard]] std::size_t n_col() const noexcept;

    // Concentrations
    [[nodiscard]] concentration_view_t get_concentration() noexcept;
    [[nodiscard]] MC::KernelConcentrationType get_device_concentration() const;

    [[nodiscard]] std::span<concentration_float_type> getConcentrationData();

    [[nodiscard]] std::span<const concentration_float_type>
    getConcentrationData() const;

    [[nodiscard]] auto getConcentrationArray() const;

    //

    [[nodiscard]] const KokkosEigen::Alias::DiagonalType<
        mass_balance_float_type>&
    getVolume() const noexcept;

    [[deprecated]] [[nodiscard]] std::span<double const> volume_span() const;

    double volume_at(std::size_t) const;

    [[nodiscard]] std::span<const double> contribution_span() const;

    [[nodiscard]] std::span<double> contribution_span_mut();

    [[nodiscard]] MC::kernelContribution get_kernel_contribution() const;

    // Setters

    void set_mass();
    // void set_kernel_contribs_to_host();
    void set_feed(std::uint64_t i_r, std::uint64_t i_c, double val);
    void set_sink(std::uint64_t i_compartment, double val);
    void set_zero_contribs();
    void setVolumes(std::span<const double> volumes,
                    std::span<const double> inv_volumes);

    void clearNegs();

  private:
    std::size_t n_r;
    std::size_t n_c;

    KokkosEigen::Alias::ColMajorMatrixtype<mass_balance_float_type> total_mass;

    KokkosEigen::KokkosEigen2D<mass_balance_float_type,
                               Kokkos::LayoutLeft,
                               Kokkos::DefaultExecutionSpace>
        concentrations;

    // FIXME: Sources does not need to available on device because this is
    // contribs which is scattered
    KokkosEigen::KokkosEigen2D<mass_balance_float_type,
                               Kokkos::LayoutRight,
                               Kokkos::DefaultExecutionSpace>
        sources;

    MC::kernelContribution contribs;
    FlowMatrixType<mass_balance_float_type> m_transition;

    KokkosEigen::Alias::DiagonalType<mass_balance_float_type> m_volumes;
    KokkosEigen::Alias::DiagonalType<mass_balance_float_type> volumes_inverse;
    KokkosEigen::Alias::DiagonalType<mass_balance_float_type> sink;

    // RowMajorEigenKokkos<double> sources;
    //

    // Kokkos::View<float**, ComputeSpace::array_layout, ComputeSpace> kc;
  };

  inline double
  ScalarSimulation::volume_at(std::size_t index) const
  {
    KOKKOS_ASSERT(index < n_c);
    return m_volumes.diagonal().coeff(EIGEN_INDEX(index));
  }

  inline auto
  ScalarSimulation::getConcentrationArray() const
  {
    // return alloc_concentrations.array();
    return concentrations.as_array();
  }

  inline MC::kernelContribution
  ScalarSimulation::get_kernel_contribution() const
  {
    // return sources.compute;
    return contribs;
  }

  inline void
  ScalarSimulation::set_zero_contribs()
  {

    sources.eigen().setZero();
    this->sink.setZero();
    sources.host_to_device_sync();
    Kokkos::deep_copy(contribs, 0);
  }

  // inline void
  // ScalarSimulation::set_kernel_contribs_to_host()
  // {
  //   // sources.update_compute_to_host();
  //   sources.device_to_host_sync();
  // }

  inline void
  ScalarSimulation::set_feed(uint64_t i_r, uint64_t i_c, double val)
  {
    this->sources.eigen().coeffRef(EIGEN_INDEX(i_r), EIGEN_INDEX(i_c)) += val;
  }

  inline void
  ScalarSimulation::set_sink(uint64_t i_compartment, double val)
  {
    this->sink.diagonal().coeffRef(EIGEN_INDEX(i_compartment)) += val;
  }

  inline const KokkosEigen::Alias::DiagonalType<mass_balance_float_type>&
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

  inline std::span<const double>
  ScalarSimulation::contribution_span() const
  {
    return this->sources.get_span();
  }

  inline std::span<double>
  ScalarSimulation::contribution_span_mut()
  {
    return this->sources.get_span();
  }

  [[deprecated]] inline std::span<double const>
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
