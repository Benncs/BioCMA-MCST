#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <cma_utils/alias.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <eigen_kokkos.hpp>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <scalar_init.hpp>
#include <scalar_simulation.hpp>
#include <simulation/descriptors/dimensions.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/kernels/move_kernel.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/probe.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <utility>
namespace Simulation
{

  // IMPORTANT to have non default to  BIND accesor(this)
  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept
      : accesor(this), mc_unit(std::move(other.mc_unit)),
        contribs_scatter(std::move(other.contribs_scatter)),
        probes(std::move(other.probes)), dims(std::move(other.dims)),
        feed(std::move(other.feed)),
        const_number_simulation(other.const_number_simulation),
        is_two_phase_flow(other.is_two_phase_flow),
        starting_time(other.starting_time), end_time(other.end_time),
        f_reaction(other.f_reaction),
        liquid_scalar(std::move(other.liquid_scalar)),
        gas_scalar(std::move(other.gas_scalar)),
        mt_model(std::move(other.mt_model)), logger(std::move(other.logger))

  {
  }

  // IMPORTANT to have non default to  BIND accesor(this)
  SimulationUnit::SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                                 ScalarInitializer&& scalar_init,
                                 std::optional<Feed::SimulationFeed> _feed)
      : accesor(this), mc_unit(std::move(_unit)),
        feed(_feed.value_or(Feed::SimulationFeed::empty())),
        is_two_phase_flow(scalar_init.gas_flow), starting_time(0), end_time(0)
  {

    this->liquid_scalar = std::make_shared<ScalarSimulation>(
        mc_unit->domain.getNumberCompartments(),
        scalar_init.n_species,
        scalar_init.volumesliq);

    this->gas_scalar = (is_two_phase_flow)
                           ? std::make_shared<ScalarSimulation>(
                                 mc_unit->domain.getNumberCompartments(),
                                 scalar_init.n_species,
                                 scalar_init.volumesgas)
                           : nullptr; // No contribs for gas

    post_init_concentration(std::move(scalar_init));

    //  scalar_init is totally moved-from here at the very end of impl
    //  post_init_concentration

    const std::size_t n_flows = this->feed.n_liquid_flow();

    mc_unit->domain.init_inner(n_flows);

    contribs_scatter
        = Kokkos::Experimental::create_scatter_view(get_kernel_contribution());

    dims = { .n_species = this->liquid_scalar->n_row(),
             .n_compartment = this->liquid_scalar->n_col() };
  }

  SimulationUnit::~SimulationUnit() = default;

  void
  SimulationUnit::reset()
  {
    liquid_scalar.reset();
    gas_scalar.reset();
    mc_unit.reset();
    // mc_unit->domain.inner = MC::DomainState<ComputeSpace, false>();
  }

  void
  SimulationUnit::updateHydro(const CmaUtils::IterationStatePtrType& newstate)
  {
    PROFILE_SECTION("simulation::updateHydro")

    const auto& liquid = newstate->get_liquid();
    const auto vl = liquid->volume();

    updateMCHydro(vl,
                  newstate->flat_neighobrs(),
                  newstate->flat_probability_leaving(),
                  liquid->out_flows());

    updateScalarHydro(newstate);
  }

  void
  SimulationUnit::updateMCHydro(std::span<const double> newliquid_volume,
                                std::span<const std::size_t> neighors_flat,
                                std::span<const double> proba_flat,
                                std::span<const double> out_flows) const
  {
    PROFILE_SECTION("simulation::updateMCHydro")
    this->mc_unit->domain.update(
        newliquid_volume, neighors_flat, out_flows, proba_flat);
  }

  void
  SimulationUnit::updateScalarHydro(
      const CmaUtils::IterationStatePtrType& newstate)
  {
    PROFILE_SECTION("simulation::updateScalarHydro")
    if (liquid_scalar) [[likely]]
    {
      const auto& liq = newstate->get_liquid();
      std::span<double const> vl = liq->volume();
      this->liquid_scalar->setVolumes(vl, liq->inverse_volume());
      this->liquid_scalar->set_transition(newstate->get_liquid()->transition());
    }
    if (newstate->has_gas() && gas_scalar)
    {
      const auto& gas = newstate->get_gas();
      this->gas_scalar->setVolumes(gas->volume(), gas->inverse_volume());
      gas_scalar->set_transition(newstate->get_gas()->transition());
      mt_model.update(newstate);
    }
  }

  void
  SimulationUnit::scatter_contribute()
  {
    auto contribs = this->get_kernel_contribution();
    Kokkos::Experimental::contribute(contribs, contribs_scatter);
    // Warning syncrho deepcopy into contribs
    //  Ok to do it here because scatter_contribute is called after
    //  cycle_process when contribs is empty (ode step is done before)
    this->liquid_scalar->synchro_sources();
  }

  void
  SimulationUnit::post_init_concentration(ScalarInitializer&& scalar_init)
  {

    if (scalar_init.type == ScalarInitialiserType::Uniform
        || scalar_init.type == ScalarInitialiserType::Local)

    {

      impl::post_init_concentration_functor(
          is_two_phase_flow, std::move(scalar_init), liquid_scalar, gas_scalar);
    }
    else
    {
      bool change_layout = scalar_init.type != ScalarInitialiserType::Serde;

      impl::post_init_concentration_file(is_two_phase_flow,
                                         std::move(scalar_init),
                                         liquid_scalar,
                                         gas_scalar,
                                         change_layout);
    }

    // At this step scalar init is totally consumed and moved-from

    // Check concentration >=0

    if ((this->liquid_scalar->getConcentrationArray() >= 0.).all())
    {
      this->liquid_scalar->set_mass();
    }
    else
    {
      throw SimulationException(ErrorCodes::BadConcentrationInitLiq);
    }

    if (is_two_phase_flow)
    {
      if ((this->gas_scalar->getConcentrationArray() >= 0.).all())
      {
        this->gas_scalar->set_mass();
      }
      else
      {
        throw SimulationException(ErrorCodes::BadConcentrationInitGas);
      }
    }
  }

} // namespace Simulation
