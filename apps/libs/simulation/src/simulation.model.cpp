#include "Kokkos_Assert.hpp"
#include <common/eigen_diag.hpp>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
EIGEN_DIAG_POP

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

  namespace
  {
    void
    set_scalar_feed(ScalarSimulation& scl, const Feed::FeedDescriptor& fd)
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

  }

  // TODO check if really noexcept
  void
  SimulationUnit::update_feed(const double d_t,
                              const bool update_scalar) noexcept
  {
    PROFILE_SECTION("host:update_feed")
    const auto rel_t = this->accesor.relative_time();
    const auto abs_t = this->accesor.absolute_time();

    auto update_feed_scalar =
        [rel_t, abs_t, d_t, update_scalar](auto& scl, Feed::FeedDescriptor& fd)
    {
      const auto time = fd.use_relative_time ? rel_t : abs_t;
      KOKKOS_ASSERT(time >= 0);
      fd.update(time, d_t);
      if (update_scalar)
      {
        set_scalar_feed(scl, fd);
      }
    };

    auto& liquid_scalar = *this->liquid_scalar;
    std::size_t mc_flow_counter = 0;

    std::ranges::for_each(
        m_feed.liquid_feeds(),
        [&](auto& feed)
        {
          update_feed_scalar(liquid_scalar, feed);

          // MC update
          if (feed.output_position.has_value())
          {
            const std::size_t output_position = *feed.output_position;
            const auto volume = liquid_scalar.volume_at(output_position);
            KOKKOS_ASSERT(volume > 0.);
            this->mc_unit->domain.set_leaving_flow(
                mc_flow_counter, output_position, feed.flow, volume);
            mc_flow_counter++;
          }
        });

    if (is_two_phase_flow)
    {
      auto& gs = *this->gas_scalar;
      auto range_feed_gas = m_feed.gas_feeds();

      std::ranges::for_each(range_feed_gas,
                            [&update_feed_scalar, &gs](auto& feed)
                            { update_feed_scalar(gs, feed); });
    }
  }

  double
  SimulationUnit::advance(double d_t) noexcept
  {
    this->m_times.advance(d_t);
    return m_times.absolute();
  }

  void
  SimulationUnit::ode_step(double d_t) const
  {

    if (is_two_phase_flow)
    {
      mt_model.gas_liquid_mass_transfer();
      const auto& mtr = mt_model.proxy()->mtr;

      this->gas_scalar->performStepGL(
          d_t, mtr, MassTransfer::Sign::GasToLiquid);

      this->liquid_scalar->performStepGL(
          d_t, mtr, MassTransfer::Sign::LiquidToGas);

      // Uncoment this to avoid negative concentration especially when low
      // oxygen concentration
      this->liquid_scalar->clearNegs();
    }
    else
    {
      this->liquid_scalar->performStep(d_t);
    }
  }
} // namespace Simulation
