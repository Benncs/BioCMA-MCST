#ifndef NDEBUG
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#  pragma GCC diagnostic ignored "-Wnan-infinity-disabled"
#endif
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifndef NDEBUG
#  pragma GCC diagnostic pop
#endif
#include <common/common.hpp>
#include <hydro/impl_mass_transfer.hpp>
#include <mc/domain.hpp>
#include <optional>
#include <scalar_simulation.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/simulation.hpp>

namespace Simulation
{

  bool
  SimulationUnit::checkScalar() const
  {
    // Ensure there is no negative values
    auto pred = [](auto&& val) { return val < 0.; };

    auto cliq = this->accesor.getCliqData();

    const auto it = std::ranges::find_if(cliq.begin(), cliq.end(), pred);
    bool flag = it == cliq.end();

    auto cgas = this->accesor.getCgasData();

    if (cgas)
    {
      const auto it = std::find_if(cgas->begin(), cgas->end(), pred);
      flag = it == cgas->end();
    }
    return flag;
  }

  // TODO check if really noexcept
  void
  SimulationUnit::clearContribution() const noexcept
  {
    if (is_two_phase_flow)
    {
      this->gas_scalar->set_zero_contribs();
    }
    this->liquid_scalar->set_zero_contribs();
  }

  // TODO check if really noexcept
  void
  SimulationUnit::update_feed(const double t,
                              const double d_t,
                              const bool update_scalar) noexcept
  {
    PROFILE_SECTION("host:update_feed")
    // Get references to the index_leaving_flow and leaving_flow data members

    auto set_scalar_feed
        = [update_scalar](auto& scl, const Feed::FeedDescriptor& fd) -> void
    {
      if (update_scalar)
      {
        // Set for each scalar, i.e for each line of the sytem, F=C_feed*Q
        for (auto [concentration, index] : fd.values)
        {
          scl.set_feed(index, fd.input_position, fd.flow * concentration);
        }
        if (fd.output_position)
        {
          scl.set_sink(*fd.output_position, fd.flow);
        }
      }
    };

    auto& liquid_scalar = *this->liquid_scalar;
    std::size_t mc_flow_counter = 0;
    for (auto& feed : feed.liquid_feeds())
    {
      feed.update(t, d_t);
      set_scalar_feed(liquid_scalar, feed);
      if (feed.output_position)
      {

        // this->mc_unit->domain.set_leaving_flow(
        //     0, *feed.output_position, feed.flow);

        this->mc_unit->domain.set_leaving_flow(
            mc_flow_counter, *feed.output_position, feed.flow);
        mc_flow_counter++;
      }
    }

    if (is_two_phase_flow)
    {
      auto& gs = *this->gas_scalar;
      for (auto& feed : feed.gas_feeds())
      {
        feed.update(t, d_t);
        set_scalar_feed(gs, feed);
      }
    }
  }

  void
  SimulationUnit::step(double d_t) const
  {

    if (is_two_phase_flow)
    {
      mt_model.gas_liquid_mass_transfer();
      const auto& mtr = mt_model.proxy()->mtr;

      this->gas_scalar->performStepGL(
          d_t, mtr, MassTransfer::Sign::GasToLiquid);

      this->liquid_scalar->performStepGL(
          d_t, mtr, MassTransfer::Sign::LiquidToGas);
    }
    else
    {
      this->liquid_scalar->performStep(d_t);
    }
  }
} // namespace Simulation
