#include <Kokkos_Core.hpp>
#include <cma_utils/iteration_state.hpp>
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
#include <scalar_simulation.hpp>
#include <simulation/alias.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/probe.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <simulation/simulation_kernel.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
#include <utility>

namespace Simulation
{

  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept = default;

  SimulationUnit::SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                                 const ScalarInitializer& scalar_init,
                                 std::optional<Feed::SimulationFeed> _feed)
      : mc_unit(std::move(_unit)), internal_counter_dead(0),
        feed(_feed.value_or(Feed::SimulationFeed{std::nullopt, std::nullopt})),
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

    Kokkos::resize(move_info.index_leaving_flow, n_flows);
    Kokkos::resize(move_info.leaving_flow, n_flows);
    Kokkos::resize(move_info.liquid_volume,
                   mc_unit->domain.getNumberCompartments());
    Kokkos::resize(move_info.diag_transition,
                   mc_unit->domain.getNumberCompartments());

           contribs_scatter=         MC::ContributionView (get_kernel_contribution());
  }

  void SimulationUnit::update(CmaUtils::IterationState&& newstate)
  {
    state = std::move(newstate);

    setVolumes();
    mc_unit->domain.setLiquidNeighbors(state.neighbors);
  }

  void SimulationUnit::scatter_contribute()
  {
    auto contribs = get_kernel_contribution();
    Kokkos::Experimental::contribute(contribs, contribs_scatter);
  }

  void SimulationUnit::setVolumes() const
  {

    std::span<double const> vl = state.liq->volume;
    std::span<double const> vg = state.gas->volume;
    if (liquid_scalar)
    {
      std::span<double const> vl = liquid_scalar->getVolumeData();
      this->liquid_scalar->setVolumes(vl, state.liq->inverse_volume);
    }

    if (gas_scalar)
    {
      this->gas_scalar->setVolumes(vg, state.gas->inverse_volume);
    }
    else
    {
      vg = std::vector<double>(vl.size(), 0);
    }

    Kokkos::View<double*, HostSpace> hostli(state.liq->volume.data(),
                                            vl.size());

    Kokkos::deep_copy(this->move_info.liquid_volume, hostli);

    this->mc_unit->domain.setVolumes(vl);
  }

  void SimulationUnit::post_init_concentration_file(
      const ScalarInitializer& scalar_init)
  {

    if (!scalar_init.liquid_buffer.has_value())
    {
      throw SimulationException(ErrorCodes::BadInitialiser);
    }

    if (!this->liquid_scalar->deep_copy_concentration(
            *scalar_init.liquid_buffer))
    {

      throw SimulationException(ErrorCodes::MismatchSize);
    }

    if (is_two_phase_flow)
    {
      if (!scalar_init.gas_buffer.has_value())
      {
        throw SimulationException(ErrorCodes::BadInitialiser);
      }
      if (!this->gas_scalar->deep_copy_concentration(*scalar_init.gas_buffer))
      {
        throw SimulationException(ErrorCodes::MismatchSize);
      }
    }
  }

  void SimulationUnit::clear_mc()
  {
    mc_unit.reset();
  }

  void SimulationUnit::reset()
  {
    liquid_scalar.reset();
    gas_scalar.reset();
    move_info = KernelInline::MoveInfo<ComputeSpace>();
  }

  void SimulationUnit::post_init_concentration_functor(
      const ScalarInitializer& scalar_init)
  {

    auto& cliq = this->liquid_scalar->get_concentration();
    decltype(&cliq) cgas = nullptr; // FIXME

    if (is_two_phase_flow)
    {
      assert(this->gas_scalar != nullptr);
      cgas = &this->gas_scalar->get_concentration();
      assert(cgas->size() != 0);
    }

    const auto& fv = scalar_init.liquid_f_init.value();

    for (decltype(cliq.rows()) i_row = 0; i_row < cliq.rows(); ++i_row)
    {
      for (decltype(cliq.cols()) i_col = 0; i_col < cliq.cols(); ++i_col)
      {
        cliq(i_row, i_col) = fv(i_row, i_col);
        if (is_two_phase_flow)
        {
          (*cgas)(i_row, i_col) = scalar_init.gas_f_init.value()(i_row, i_col);
        }
      }
    }
  }

  void
  SimulationUnit::post_init_concentration(const ScalarInitializer& scalar_init)
  {

    if (scalar_init.type == ScalarInitialiserType::Uniform ||
        scalar_init.type == ScalarInitialiserType::Local)

    {
      post_init_concentration_functor(scalar_init);
    }
    else
    {
      post_init_concentration_file(scalar_init);
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
