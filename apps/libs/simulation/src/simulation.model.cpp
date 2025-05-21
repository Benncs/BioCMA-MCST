#include "simulation/mass_transfer.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <common/common.hpp>
#include <hydro/impl_mass_transfer.hpp>
#include <mc/domain.hpp>
#include <optional>
#include <scalar_simulation.hpp>
#include <simulation/simulation.hpp>

namespace Simulation
{

  std::span<const double> SimulationUnit::getContributionData() const
  {
    return liquid_scalar->getContributionData();
  }

  std::span<double> SimulationUnit::getContributionData_mut()
  {
    return liquid_scalar->getContributionData();
  }

  std::span<double> SimulationUnit::getCliqData() const
  {
    return this->liquid_scalar->getConcentrationData();
  }

  [[nodiscard]] std::optional<std::span<const double>> SimulationUnit::getCgasData() const
  {
    if (!gas_scalar)
    {
      return std::nullopt;
    }
    return this->gas_scalar->getConcentrationData();
  }

  [[nodiscard]] Dimensions SimulationUnit::getDimensions() const noexcept
  {
    return {this->liquid_scalar->n_row(), this->liquid_scalar->n_col()};
  }

  [[nodiscard]] std::optional<std::span<const double>> SimulationUnit::getMTRData() const
  {
    return this->mt_model.mtr_data();
  }

  [[deprecated("perf:not useful")]] void
  SimulationUnit::reduceContribs_per_rank(std::span<const double> data) const
  {

    PROFILE_SECTION("host:reduceContribs_rank")
    this->liquid_scalar->reduce_contribs(data);
  }

  void SimulationUnit::reduceContribs(std::span<const double> data, size_t n_rank) const
  {
    PROFILE_SECTION("host:reduceContribs")
    const auto [nr, nc] = getDimensions();
    this->liquid_scalar->set_zero_contribs();

    for (int i = 0; i < static_cast<int>(n_rank); ++i)
    {
      this->liquid_scalar->reduce_contribs({&data[i * nr * nc], nr * nc});
    }
  }

  void SimulationUnit::clearContribution() const noexcept
  {
    // this->liquid_scalar->vec_kla.setZero();
    // Dont forget to clear kernel contribution
    if (is_two_phase_flow)
    {
      this->gas_scalar->set_zero_contribs();
    }
    this->liquid_scalar->set_zero_contribs();
  }

  void
  SimulationUnit::update_feed(const double t, const double d_t, const bool update_scalar) noexcept
  {
    PROFILE_SECTION("host:update_feed")
    // Get references to the index_leaving_flow and leaving_flow data members
    const auto& _index_leaving_flow = this->move_info.index_leaving_flow;
    const auto& _leaving_flow = this->move_info.leaving_flow;

    // Get the index of the exit compartment
    // TODO exit is not necessarly at the index n-1, it should be given by user
    const uint64_t i_exit = 0; // mc_unit->domain.getNumberCompartments() - 1;

    // Define the set_feed lambda function
    auto set_feed =
        [t, d_t, i_exit, &_index_leaving_flow, &_leaving_flow, update_scalar](
            const std::shared_ptr<ScalarSimulation>& scalar, auto&& descritor, bool mc_f = false)

    {
      double flow = 0.; // Initialize the flow variable
      bool set_exit = false;
      // Iterate through each current_feed in the descriptor
      for (auto&& current_feed : descritor)
      {
        current_feed.update(t, d_t);     // Update the current_feed
        flow += current_feed.flow_value; // Get the flow_value of the current_feed
        set_exit = current_feed.set_exit;
        if (update_scalar)
        {
          // Iterate through the species, positions, and values of the
          // current_feed
          for (std::size_t i_f = 0; i_f < current_feed.n_v; ++i_f)
          {
            const std::size_t i_species = current_feed.species[i_f];
            scalar->set_feed(i_species, current_feed.position[i_f], flow * current_feed.value[i_f]);
          }
        }
      }

      if (set_exit)
      {
        if (update_scalar)
        {
          // Set the sink for the exit compartment
          scalar->set_sink(i_exit, flow);
        }

        // Update Flow for mc particle
        if (mc_f)
        {
          _index_leaving_flow(0) = i_exit;
          _leaving_flow(0) = flow;
        }
      }
    };

    if (feed.liquid.has_value())
    {
      set_feed(this->liquid_scalar, *feed.liquid, true);
    }

    if (is_two_phase_flow && feed.gas.has_value())
    {
      set_feed(this->gas_scalar, *feed.gas);
    }
  }

  void SimulationUnit::step(double d_t) const
  {

    if (is_two_phase_flow)
    {
      mt_model.gas_liquid_mass_transfer(state);
      const auto& mtr = mt_model.proxy()->mtr;

      this->gas_scalar->performStepGL(
          d_t, state.gas->get_transition(), mtr, MassTransfer::Sign::GasToLiquid);

      this->liquid_scalar->performStepGL(
          d_t, state.liq->get_transition(), mtr, MassTransfer::Sign::LiquidToGas);
    }
    else
    {
      this->liquid_scalar->performStep(d_t, state.liq->get_transition());
    }
  }
} // namespace Simulation
