#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <algorithm>
#include <cma_utils/alias.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <scalar.hpp>
#include <scalar_init.hpp>
#include <simulation/descriptors/dimensions.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/kernels/move_kernel.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/probe.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <stdexcept>
#include <utility>
namespace Simulation
{

  // IMPORTANT to have non default to  BIND accesor(this)
  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept
      : accesor(this), mc_unit(std::move(other.mc_unit)),
        m_table(std::move(other.m_table)),
        contribs_scatter(std::move(other.contribs_scatter)),
        probes(std::move(other.probes)), dims(std::move(other.dims)),
        m_feed(std::move(other.m_feed)),
        const_number_simulation(other.const_number_simulation),
        is_two_phase_flow(other.is_two_phase_flow),
        f_reaction(other.f_reaction), m_times(other.m_times),
        scalar(std::move(other.scalar)),
        m_mt_model(std::move(other.m_mt_model)), logger(std::move(other.logger))

  {
  }

  // IMPORTANT to have non default to  BIND accesor(this)
  SimulationUnit::SimulationUnit(
      std::unique_ptr<MC::MonteCarloUnit>&& _unit,
      ScalarInitializer&& scalar_init,
      std::optional<Feed::SimulationFeed> _feed,
      std::shared_ptr<Mixture::SpecieTable> table,
      std::optional<MassTransfer::Type::MtrTypeVariant>&& mt_model_variant,
      Sparam params)
      : accesor(this), mc_unit(std::move(_unit)), m_table(std::move(table)),
        m_feed(_feed.value_or(Feed::SimulationFeed::empty())),
        is_two_phase_flow(scalar_init.gas_flow), f_reaction(params.f_reaction)

  {

    const auto n_compartments = mc_unit->domain.getNumberCompartments();

    this->scalar
        = (is_two_phase_flow)
              ? std::make_shared<impl::ScalarSimulation>(n_compartments,
                                                         scalar_init.n_species,
                                                         scalar_init.volumesliq,
                                                         scalar_init.volumesgas)
              : std::make_shared<impl::ScalarSimulation>(
                    n_compartments,
                    scalar_init.n_species,
                    scalar_init.volumesliq);

    post_init_concentration(std::move(scalar_init));

    if (mt_model_variant)
    {
      setMtrModel(std::move(*mt_model_variant));
    }
    //  scalar_init is totally moved-from here at the very end of impl
    //  post_init_concentration

    const std::size_t n_flows = m_feed.n_liquid_flow();

    mc_unit->domain.init_inner(n_flows);

    contribs_scatter
        = Kokkos::Experimental::create_scatter_view(get_kernel_contribution());

    dims = { .n_species = this->scalar->n_species(),
             .n_compartment = this->scalar->n_compartments() };
  }

  SimulationUnit::~SimulationUnit() = default;

  void
  SimulationUnit::reset()
  {
    scalar.reset();
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
    if (!scalar) [[unlikely]]
    {
      throw std::runtime_error("No scalar defined");
    }

    const auto& liq = newstate->get_liquid();
    this->scalar->set_volumes(
        Phase::Liquid, liq->volume(), liq->inverse_volume());

    if (scalar->has_gas())
    {
      const auto& gas = newstate->get_gas();
      this->scalar->set_volumes(
          Phase::Gas, gas->volume(), gas->inverse_volume());
      this->scalar->set_transition_two_phase(liq->transition(),
                                             gas->transition());

      m_mt_model.update(newstate);
    }
    else
    {
      this->scalar->set_transition(liq->transition());
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
    this->scalar->synchro_sources();
  }

  void
  SimulationUnit::post_init_concentration(ScalarInitializer&& scalar_init)
  {

    if (scalar_init.type == ScalarInitialiserType::Uniform
        || scalar_init.type == ScalarInitialiserType::Local)

    {

      impl::post_init_concentration_functor(
          is_two_phase_flow, std::move(scalar_init), scalar);
    }
    else
    {
      bool change_layout = scalar_init.type != ScalarInitialiserType::Serde;

      impl::post_init_concentration_file(
          is_two_phase_flow, std::move(scalar_init), scalar, change_layout);
    }

    // At this step scalar init is totally consumed and moved-from

    // Check concentration >=0

    auto positive = [](double c) { return c >= 0.; };

    if (!std::ranges::all_of(this->scalar->concentration_data(Phase::Liquid),
                             positive))
    {
      throw BaseSimulationException(ErrorCodes::BadConcentrationInitLiq);
    }

    if (is_two_phase_flow
        && !std::ranges::all_of(this->scalar->concentration_data(Phase::Gas),
                                positive))
    {
      throw BaseSimulationException(ErrorCodes::BadConcentrationInitGas);
    }

    this->scalar->set_mass();
  }

} // namespace Simulation
