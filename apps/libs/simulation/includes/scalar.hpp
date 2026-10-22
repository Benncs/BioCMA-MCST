#ifndef __SCALAR_SIMULATION_TP_HPP__
#define __SCALAR_SIMULATION_TP_HPP__

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
#include <optional>
#include <simulation/mass_transfer.hpp>
#include <span>

namespace Simulation
{
  using mass_balance_float_type = double;

  namespace impl
  {
    using FlowMatrix = Eigen::SparseMatrix<mass_balance_float_type>;

    using Mat = KokkosEigen::Alias::ColMajorMatrixtype<mass_balance_float_type>;

    using Diag = KokkosEigen::Alias::DiagonalType<mass_balance_float_type>;

    using KokkosMat = KokkosEigen::KokkosEigen2D<mass_balance_float_type,
                                                 Kokkos::LayoutLeft,
                                                 Kokkos::DefaultExecutionSpace>;

    // Column blocks of ColMajor storage are contiguous (InnerPanel=true)
    using MatBlock = Mat::ColsBlockXpr;
    using MapBlock = KokkosMat::eigen_map_type::ColsBlockXpr;
    using VecBlock = std::remove_reference_t<
        decltype(std::declval<Diag&>().diagonal())>::SegmentReturnType;

    /**
     * @brief Augmented system: phases are stacked by columns,
     * every matrix is n_species x (n_phase * n_compartments) = [ L | G ].
     *
     * With transition = blockdiag(M_L, M_G) the whole explicit step is a
     * single expression.
     */
    struct SystemScalar
    {
      Mat total_mass;
      KokkosMat concentrations;

      KokkosMat sources;
      // [ mtr | -mtr ], stays 0 for one phase
      Mat mtr;

      FlowMatrix transition;

      Diag volumes;
      Diag volumes_inverse;
      Diag sink;

      SystemScalar(std::size_t n_species, std::size_t n_cols);

      SystemScalar(const SystemScalar&) = delete;
      SystemScalar(SystemScalar&&) = delete;
      SystemScalar& operator=(const SystemScalar&) = delete;
      SystemScalar& operator=(SystemScalar&&) = delete;
      ~SystemScalar() = default;

      void step_explicit(double d_t);

      void update_concentration();
    };

    /**
     * @brief Views on one phase block of SystemScalar, built once.
     *
     * Valid as long as SystemScalar does not move (non movable) and is never
     * resized (by construction).
     */
    struct PhaseView
    {
      Phase identifier;
      MapBlock concentration;
      MatBlock mass;
      MapBlock source;
      MatBlock mtr;
      VecBlock volume;
      VecBlock volume_inverse;
      VecBlock sink;

      PhaseView(Phase id,
                SystemScalar& system,
                std::size_t n_species,
                std::size_t n_compartments);

      PhaseView(const PhaseView&) = delete;
      PhaseView(PhaseView&&) = delete;
      PhaseView& operator=(const PhaseView&) = delete;
      PhaseView& operator=(PhaseView&&) = delete;
      ~PhaseView() = default;
    };

    class ScalarSimulation
    {
    public:
      using concentration_float_type = mass_balance_float_type;

      // Liquid only
      ScalarSimulation(std::size_t n_compartments,
                       std::size_t n_species,
                       std::span<const mass_balance_float_type> volume_liq);

      // Gas-liquid
      ScalarSimulation(std::size_t n_compartments,
                       std::size_t n_species,
                       std::span<const mass_balance_float_type> volume_liq,
                       std::span<const mass_balance_float_type> volume_gas);

      ScalarSimulation(const ScalarSimulation&) = delete;
      ScalarSimulation(ScalarSimulation&&) = delete;
      ScalarSimulation& operator=(const ScalarSimulation&) = delete;
      ScalarSimulation& operator=(ScalarSimulation&&) = delete;
      ~ScalarSimulation() = default;

      // Hydro

      void set_transition(CmaUtils::StateCooMatrixType&& transition_liq);

      void
      set_transition_two_phase(CmaUtils::StateCooMatrixType&& transition_liq,
                               CmaUtils::StateCooMatrixType&& transition_gas);

      void set_volumes(Phase phase,
                       std::span<const double> volumes,
                       std::span<const double> inv_volumes);

      // Step

      void step_explicit(double d_t);

      void step_explicit_two_phase(
          double d_t, const MassTransfer::MassTransferProxy& mass_transfer);

      void clear_negs();

      // Sources

      void set_feed(Phase phase,
                    std::size_t i_species,
                    std::size_t i_compartment,
                    double val);

      void set_sink(Phase phase, std::size_t i_compartment, double val);

      void set_zero_contribs();

      void synchro_sources();

      // Init

      bool deep_copy_concentration(Phase phase, std::span<const double> data);

      void set_mass();

      // Getters

      [[nodiscard]] std::span<const double>
      concentration_data(Phase phase) const;

      [[nodiscard]] std::span<double> concentration_data(Phase phase);

      // Host Kokkos view on one phase block (LayoutLeft, contiguous)
      [[nodiscard]] KokkosMat::host_view_type
      concentration_host_view(Phase phase);

      [[nodiscard]] std::span<const double> source_data(Phase phase) const;

      [[nodiscard]] std::span<double> source_data(Phase phase);

      [[nodiscard]] double volume_at(Phase phase,
                                     std::size_t i_compartment) const;

      [[nodiscard]] MC::KernelConcentrationType
      get_device_concentration() const;

      [[nodiscard]] MC::kernelContribution get_kernel_contribution() const;

      [[nodiscard]] std::size_t n_species() const noexcept;

      [[nodiscard]] std::size_t n_compartments() const noexcept;

      [[nodiscard]] bool has_gas() const noexcept;

      // Views given once to MassTransferModel
      [[nodiscard]] MassTransfer::MtrInput mtr_input();

    private:
      [[nodiscard]] PhaseView& phase(Phase phase);
      [[nodiscard]] const PhaseView& phase(Phase phase) const;
      [[nodiscard]] bool is_valid_phase(Phase phase) const;

      [[nodiscard]] std::size_t offset(Phase phase) const noexcept;

      std::size_t n_s;
      std::size_t n_c;
      std::size_t n_phase;

      SystemScalar system;
      PhaseView liquid;
      std::optional<PhaseView> gas;

      MC::kernelContribution contribs;
    };

  } // namespace impl

} // namespace Simulation

#endif //__SCALAR_SIMULATION_TP_HPP__
