#include "cma_utils/d_transitionner.hpp"
#include "common/common.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
// #include <cma_utils/iteration_state.hpp>
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
#include <simulation/alias.hpp>
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

  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept = default;

  SimulationUnit::SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                                 const ScalarInitializer& scalar_init,
                                 std::optional<Feed::SimulationFeed> _feed)
      : mc_unit(std::move(_unit)),
        feed(_feed.value_or(Feed::SimulationFeed::empty())),
        is_two_phase_flow(scalar_init.gas_flow)
  {

    this->liquid_scalar = std::make_shared<ScalarSimulation>(
        mc_unit->domain.getNumberCompartments(),
        scalar_init.n_species,
        scalar_init.volumesliq);

    this->gas_scalar = (is_two_phase_flow)
                           ? std::make_shared<ScalarSimulation>(
                                 mc_unit->domain.getNumberCompartments(),
                                 scalar_init.n_species,
                                 scalar_init.volumesgas) // No contribs for gas
                           : nullptr;

    post_init_concentration(scalar_init);
    post_init_compartments();

    const std::size_t n_flows = this->feed.n_liquid_flow();

    move_info = KernelInline::MoveInfo<ComputeSpace>(
        mc_unit->domain.getNumberCompartments(), n_flows);

    contribs_scatter =
        Kokkos::Experimental::create_scatter_view(get_kernel_contribution());
  }

  // void SimulationUnit::update(CmaUtils::IterationState&& new_state)
  // {
  //   state = std::move(new_state);
  //   setVolumes();
  //   mc_unit->domain.setLiquidNeighbors(state.neighbors);
  // }

  void
  SimulationUnit::updateHydro(const CmaUtils::IterationStatePtrType& newstate)
  {
    PROFILE_SECTION("simulation::updateHydro")

    auto liquid = newstate->get_liquid();
    auto vl = liquid->volume();

    updateMCHydro(vl,
                  newstate->flat_neighobrs(),
                  newstate->flat_probability_leaving(),
                  liquid->out_flows());

    updateScalarHydro(newstate);
  }

  void SimulationUnit::updateMCHydro(std::span<const double> vl,
                                     std::span<const std::size_t> neighors_flat,
                                     std::span<const double> proba_flat,
                                     std::span<const double> out_flows)
  {
    PROFILE_SECTION("simulation::updateMCHydro")
    this->mc_unit->domain.setVolumes(vl);
    this->mc_unit->domain.setLiquidNeighbors(neighors_flat);

    this->move_info.update(vl, out_flows, proba_flat);
  }

  void SimulationUnit::updateScalarHydro(
      const CmaUtils::IterationStatePtrType& newstate)
  {
    PROFILE_SECTION("simulation::updateScalarHydro")
    this->setVolumes(newstate);

    if (liquid_scalar)
    {
      liquid_scalar->set_transition(newstate->get_liquid()->transition());
    }
    if (newstate->has_gas() && gas_scalar)
    {
      gas_scalar->set_transition(newstate->get_gas()->transition());
      mt_model.update(newstate);
    }
  }

  void SimulationUnit::scatter_contribute()
  {
    auto contribs = get_kernel_contribution();
    Kokkos::Experimental::contribute(contribs, contribs_scatter);
    // Warning syncrho deepcopy into contribs
    //  Ok to do it here because scatter_contribute is called after
    //  cycle_process when contribs is empty (ode step is done before)
    this->liquid_scalar->synchro_sources();
  }

  void SimulationUnit::setVolumes(
      const CmaUtils::IterationStatePtrType& newstate) const
  {
    auto liq = newstate->get_liquid();
    std::span<double const> vl = liq->volume();
    if (liquid_scalar)
    {
      this->liquid_scalar->setVolumes(vl, liq->inverse_volume());
    }

    if (newstate->has_gas() && gas_scalar)
    {
      auto gas = newstate->get_gas();
      this->gas_scalar->setVolumes(gas->volume(), gas->inverse_volume());
    }
  }

  void SimulationUnit::reset()
  {
    liquid_scalar.reset();
    gas_scalar.reset();
    move_info = KernelInline::MoveInfo<ComputeSpace>();
  }

  void
  SimulationUnit::post_init_concentration(const ScalarInitializer& scalar_init)
  {

    if (scalar_init.type == ScalarInitialiserType::Uniform ||
        scalar_init.type == ScalarInitialiserType::Local)

    {

      impl::post_init_concentration_functor(
          is_two_phase_flow, scalar_init, liquid_scalar, gas_scalar);
    }
    else
    {
      impl::post_init_concentration_file(
          is_two_phase_flow, scalar_init, liquid_scalar, gas_scalar);
    }

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

  void SimulationUnit::post_init_compartments()
  {
    // auto _compute_concentration = liquid_scalar->get_device_concentration();
    // // auto _containers = mc_unit->domain.data();

    // Kokkos::parallel_for(
    //     "post_init_compartments",
    //     Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,
    //                                                        mc_unit->domain.getNumberCompartments()),
    //     KOKKOS_LAMBDA(const int i) {
    //       // auto s = Kokkos::subview(_compute_concentration, i,
    //       Kokkos::ALL); _containers(i).concentrations =
    //       Kokkos::subview(_compute_concentration, Kokkos::ALL, i);
    //     });
    // Kokkos::fence();
  }

  SimulationUnit::~SimulationUnit() = default;

} // namespace Simulation
