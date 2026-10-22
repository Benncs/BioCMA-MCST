#include <Kokkos_Assert.hpp>
#include <common/eigen_diag.hpp>
#include <mc/alias.hpp>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <kokkos_eigen.hpp>
EIGEN_DIAG_POP
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <common/common.hpp>
#include <hydro/impl_mass_transfer.hpp>
#include <scalar.hpp>
#include <simulation/simulation_exception.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

#define CHECK_PHASE(__phase__)                                                 \
  if (!is_valid_phase(phase))                                                  \
  {                                                                            \
    throw std::invalid_argument("ScalarSimulation: no gas phase");             \
  }

namespace
{
  using Triplet = Eigen::Triplet<double>;

  std::size_t
  phase_column(Phase phase, std::size_t n_compartments) noexcept
  {
    return (phase == Phase::Gas) ? n_compartments : 0;
  }

  void
  append_triplets(std::vector<Triplet>& triplets,
                  const CmaUtils::StateCooMatrixType& transition,
                  std::size_t offset)
  {
    const auto rows = transition->row_indices();
    const auto cols = transition->col_indices();
    const auto vals = transition->values();
    KOKKOS_ASSERT(rows.size() == vals.size() && cols.size() == vals.size());

    for (std::size_t i = 0; i < vals.size(); ++i)
    {
      triplets.emplace_back(EIGEN_INDEX(rows[i] + offset),
                            EIGEN_INDEX(cols[i] + offset),
                            vals[i]);
    }
  }

  void
  copy_volume(Simulation::impl::VecBlock& dst, std::span<const double> src)
  {
    if (std::cmp_not_equal(src.size(), dst.size()))
    {
      throw Simulation::BaseSimulationException(
          Simulation::ErrorCodes::MismatchSizeVolume);
    }
    std::ranges::copy(src, dst.begin());
  }

} // namespace

namespace Simulation::impl
{

  // SystemScalar

  SystemScalar::SystemScalar(std::size_t n_species, std::size_t n_cols)
      : total_mass(n_species, n_cols),
        concentrations("concentrations", n_species, n_cols),
        sources("sources", n_species, n_cols), mtr(n_species, n_cols),
        transition(EIGEN_INDEX(n_cols), EIGEN_INDEX(n_cols)),
        volumes(EIGEN_INDEX(n_cols)), volumes_inverse(EIGEN_INDEX(n_cols)),
        sink(EIGEN_INDEX(n_cols))
  {
    total_mass.setZero();
    concentrations.eigen().setZero();
    sources.eigen().setZero();
    mtr.setZero();
    volumes.setZero();
    volumes_inverse.setIdentity();
    sink.setZero();
  }

  void
  SystemScalar::step_explicit(double d_t)
  {
    PROFILE_SECTION("performStep_augmented")

    const auto& c = concentrations.cst_eigen();

    total_mass.noalias()
        += d_t * (c * transition - c * sink + sources.cst_eigen() + mtr);

    update_concentration();
  }

  void
  SystemScalar::update_concentration()
  {
    concentrations.eigen().noalias() = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.host_to_device_sync();
  }

  // PhaseView

  PhaseView::PhaseView(Phase id,
                       SystemScalar& system,
                       std::size_t n_species,
                       std::size_t n_compartments)
      : identifier(id), concentration(system.concentrations.eigen().middleCols(
                            EIGEN_INDEX(phase_column(id, n_compartments)),
                            EIGEN_INDEX(n_compartments))),
        mass(system.total_mass.middleCols(
            EIGEN_INDEX(phase_column(id, n_compartments)),
            EIGEN_INDEX(n_compartments))),
        source(system.sources.eigen().middleCols(
            EIGEN_INDEX(phase_column(id, n_compartments)),
            EIGEN_INDEX(n_compartments))),
        mtr(system.mtr.middleCols(EIGEN_INDEX(phase_column(id, n_compartments)),
                                  EIGEN_INDEX(n_compartments))),
        volume(system.volumes.diagonal().segment(
            EIGEN_INDEX(phase_column(id, n_compartments)),
            EIGEN_INDEX(n_compartments))),
        volume_inverse(system.volumes_inverse.diagonal().segment(
            EIGEN_INDEX(phase_column(id, n_compartments)),
            EIGEN_INDEX(n_compartments))),
        sink(system.sink.diagonal().segment(
            EIGEN_INDEX(phase_column(id, n_compartments)),
            EIGEN_INDEX(n_compartments)))
  {
    KOKKOS_ASSERT(std::cmp_equal(concentration.rows(), n_species));
  }

  // ScalarSimulation

  ScalarSimulation::ScalarSimulation(
      std::size_t n_compartments,
      std::size_t n_species,
      std::span<const mass_balance_float_type> volume_liq)
      : n_s(n_species), n_c(n_compartments), n_phase(1),
        system(n_s, n_phase * n_c), liquid(Phase::Liquid, system, n_s, n_c),
        contribs("contribs", n_s, n_c)
  {
    copy_volume(liquid.volume, volume_liq);
  }

  ScalarSimulation::ScalarSimulation(
      std::size_t n_compartments,
      std::size_t n_species,
      std::span<const mass_balance_float_type> volume_liq,
      std::span<const mass_balance_float_type> volume_gas)
      : n_s(n_species), n_c(n_compartments), n_phase(2),
        system(n_s, n_phase * n_c), liquid(Phase::Liquid, system, n_s, n_c),
        gas(std::in_place, Phase::Gas, system, n_s, n_c),
        contribs("contribs", n_s, n_c)
  {
    copy_volume(liquid.volume, volume_liq);
    copy_volume(gas->volume, volume_gas);
  }

  // Hydro

  void
  ScalarSimulation::set_transition(
      CmaUtils::StateCooMatrixType&& transition_liq)
  {
    if (has_gas())
    {
      throw std::invalid_argument(
          "set_transition: two phase system needs set_transition_two_phase");
    }
    std::vector<Triplet> triplets;
    triplets.reserve(transition_liq->values().size());
    append_triplets(triplets, transition_liq, 0);

    system.transition.setFromTriplets(triplets.begin(), triplets.end());
    system.transition.makeCompressed();
  }

  void
  ScalarSimulation::set_transition_two_phase(
      CmaUtils::StateCooMatrixType&& transition_liq,
      CmaUtils::StateCooMatrixType&& transition_gas)
  {
    if (!has_gas())
    {
      throw std::invalid_argument(
          "set_transition_two_phase: system has no gas phase");
    }
    std::vector<Triplet> triplets;
    triplets.reserve(transition_liq->values().size()
                     + transition_gas->values().size());
    append_triplets(triplets, transition_liq, 0);
    // Gas block on the diagonal: transition = blockdiag(M_L, M_G)
    append_triplets(triplets, transition_gas, n_c);

    system.transition.setFromTriplets(triplets.begin(), triplets.end());
    system.transition.makeCompressed();
  }

  void
  ScalarSimulation::set_volumes(Phase phase,
                                std::span<const double> volumes,
                                std::span<const double> inv_volumes)
  {
    auto& view = this->phase(phase);
    copy_volume(view.volume, volumes);
    copy_volume(view.volume_inverse, inv_volumes);
  }

  // Step

  void
  ScalarSimulation::step_explicit(double d_t)
  {
    system.step_explicit(d_t);
  }

  void
  ScalarSimulation::step_explicit_two_phase(
      double d_t, const MassTransfer::MassTransferProxy& mass_transfer)
  {
    if (!has_gas())
    {
      throw std::invalid_argument(
          "step_explicit_two_phase: system has no gas phase");
    }
    // Same convention as MassTransfer::Sign: LiquidToGas=+1 on liquid,
    // GasToLiquid=-1 on gas
    liquid.mtr = mass_transfer.mtr;
    gas->mtr = -mass_transfer.mtr;

    system.step_explicit(d_t);
  }

  void
  ScalarSimulation::clear_negs()
  {
    auto& concentrations = system.concentrations;
    using float_type = KokkosMat::float_type;
    // order of magniture species to clip
    constexpr float_type max_species_value = 5e-3;
    // O(d_t)=1e-4
    constexpr float_type scheme_relative_error = 1e-4;

    constexpr float_type TOL = scheme_relative_error * max_species_value;

    using space = KokkosMat::host_view_type::execution_space;

    // Liquid only: columns [0, n_c)
    auto hv = concentrations.host_view();
    Kokkos::parallel_for(
        "clear_negs",
        Kokkos::MDRangePolicy<
            space,
            Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
            { 0, 0 }, { n_s, n_c }),
        KOKKOS_LAMBDA(int i, int j) {
          const auto val = hv(i, j);
          if (val < static_cast<float_type>(0) && Kokkos::abs(val) < TOL)
          {
            hv(i, j) = static_cast<float_type>(0);
          }
        });

    concentrations.host_to_device_sync();
  }

  // Sources

  void
  ScalarSimulation::set_feed(Phase phase,
                             std::size_t i_species,
                             std::size_t i_compartment,
                             double val)
  {
    this->phase(phase).source.coeffRef(EIGEN_INDEX(i_species),
                                       EIGEN_INDEX(i_compartment))
        += val;
  }

  void
  ScalarSimulation::set_sink(Phase phase, std::size_t i_compartment, double val)
  {
    this->phase(phase).sink.coeffRef(EIGEN_INDEX(i_compartment)) += val;
  }

  void
  ScalarSimulation::set_zero_contribs()
  {
    system.sources.eigen().setZero();
    system.sink.setZero();
    system.sources.host_to_device_sync();
    Kokkos::deep_copy(contribs, 0);
  }

  void
  ScalarSimulation::synchro_sources()
  {
    // Only liquid columns come from MC, gas columns on host are left untouched
    const auto liquid_cols = std::make_pair(std::size_t{ 0 }, n_c);
    auto device_liquid = Kokkos::subview(
        system.sources.device_view(), Kokkos::ALL, liquid_cols);
    auto host_liquid
        = Kokkos::subview(system.sources.host_view(), Kokkos::ALL, liquid_cols);

    Kokkos::deep_copy(device_liquid, contribs);
    Kokkos::deep_copy(host_liquid, device_liquid);
  }

  // Init

  bool
  ScalarSimulation::deep_copy_concentration(Phase phase,
                                            std::span<const double> data)
  {
    auto dst = concentration_data(phase);
    if (data.size() != dst.size())
    {
      return false;
    }
    std::ranges::copy(data, dst.begin());
    system.concentrations.host_to_device_sync();
    return true;
  }

  void
  ScalarSimulation::set_mass()
  {
    system.total_mass = system.concentrations.cst_eigen() * system.volumes;
  }

  // Getters

  std::span<const double>
  ScalarSimulation::concentration_data(Phase phase) const
  {
    CHECK_PHASE(phase)
    return system.concentrations.get_span().subspan(offset(phase) * n_s,
                                                    n_s * n_c);
  }

  std::span<double>
  ScalarSimulation::concentration_data(Phase phase)
  {

    CHECK_PHASE(phase)

    return system.concentrations.get_span().subspan(offset(phase) * n_s,
                                                    n_s * n_c);
  }

  KokkosMat::host_view_type
  ScalarSimulation::concentration_host_view(Phase phase)
  {
    CHECK_PHASE(phase)
    return Kokkos::subview(system.concentrations.host_view(),
                           Kokkos::ALL,
                           std::make_pair(offset(phase), offset(phase) + n_c));
  }

  std::span<const double>
  ScalarSimulation::source_data(Phase phase) const
  {
    CHECK_PHASE(phase)
    return system.sources.get_span().subspan(offset(phase) * n_s, n_s * n_c);
  }

  std::span<double>
  ScalarSimulation::source_data(Phase phase)
  {
    CHECK_PHASE(phase)
    return system.sources.get_span().subspan(offset(phase) * n_s, n_s * n_c);
  }

  double
  ScalarSimulation::volume_at(Phase phase, std::size_t i_compartment) const
  {
    KOKKOS_ASSERT(i_compartment < n_c);
    return this->phase(phase).volume.coeff(EIGEN_INDEX(i_compartment));
  }

  MC::KernelConcentrationType
  ScalarSimulation::get_device_concentration() const
  {
    // Liquid columns only: contiguous in LayoutLeft, stays LayoutLeft
    return Kokkos::subview(system.concentrations.device_view_cst(),
                           Kokkos::ALL,
                           std::make_pair(std::size_t{ 0 }, n_c));
  }

  MC::kernelContribution
  ScalarSimulation::get_kernel_contribution() const
  {
    return contribs;
  }

  std::size_t
  ScalarSimulation::n_species() const noexcept
  {
    return n_s;
  }

  std::size_t
  ScalarSimulation::n_compartments() const noexcept
  {
    return n_c;
  }

  bool
  ScalarSimulation::has_gas() const noexcept
  {
    return gas.has_value();
  }

  [[nodiscard]] bool
  ScalarSimulation::is_valid_phase(Phase phase) const
  {
    if (phase == Phase::Liquid)
    {
      return true;
    }
    return has_gas();
  }

  const PhaseView&
  ScalarSimulation::phase(Phase phase) const
  {
    if (phase == Phase::Liquid)
    {
      return liquid;
    }
    if (!has_gas())
    {
      throw std::invalid_argument("ScalarSimulation: no gas phase");
    }
    return *gas;
  }

  PhaseView&
  ScalarSimulation::phase(Phase phase)
  {
    if (phase == Phase::Liquid)
    {
      return liquid;
    }
    if (!has_gas())
    {
      throw std::invalid_argument("ScalarSimulation: no gas phase");
    }
    return *gas;
  }

  std::size_t
  ScalarSimulation::offset(Phase phase) const noexcept
  {
    return phase_column(phase, n_c);
  }

  MassTransfer::MtrInput
  ScalarSimulation::mtr_input()
  {
    MassTransfer::MtrInput input
        = { concentration_host_view(Phase::Liquid),
            concentration_host_view(Phase::Gas),
            MassTransfer::MtrVolumeView(liquid.volume.data(), n_c) };

    return input;
  }

} // namespace Simulation::impl
