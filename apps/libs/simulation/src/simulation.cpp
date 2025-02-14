#include "simulation/mass_transfer.hpp"
#include "simulation/simulation_kernel.hpp"
#include <Kokkos_Core.hpp>
#include <cma_read/light_2d_view.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <scalar_simulation.hpp>
#include <simulation/alias.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
#include <cma_utils/transport.hpp>
#include <utility>

namespace Simulation
{

  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept = default;
  //     : mc_unit(std::move(other.mc_unit)), internal_counter_dead(other.internal_counter_dead),
  //       waiting_allocation_particle(std::move(other.waiting_allocation_particle)),
  //       flow_liquid(other.flow_liquid), flow_gas(other.flow_gas),
  //       is_two_phase_flow(other.is_two_phase_flow)
  // {
  // }

  SimulationUnit::SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                                 const ScalarInitializer& scalar_init,
                                 std::optional<Feed::SimulationFeed> _feed)
      : mc_unit(std::move(_unit)), internal_counter_dead(0),
        waiting_allocation_particle("waiting_allocation_particle"),
        feed(_feed.value_or(Feed::SimulationFeed{std::nullopt, std::nullopt})),
        is_two_phase_flow(scalar_init.gas_flow), move_info()
  {

    this->liquid_scalar = std::make_shared<ScalarSimulation>(
        mc_unit->domain.getNumberCompartments(), scalar_init.n_species, scalar_init.volumesliq);

    this->gas_scalar =
        (is_two_phase_flow)
            ? std::make_shared<ScalarSimulation>(mc_unit->domain.getNumberCompartments(),
                                                 scalar_init.n_species,
                                                 scalar_init.volumesgas) // No contribs for gas
            : nullptr;

    post_init_concentration(scalar_init);
    post_init_compartments();

    const std::size_t n_flows = (this->feed.liquid.has_value()) ? this->feed.liquid->size() : 0;

    move_info.index_leaving_flow = LeavingFlowIndexType("index_leaving_flow", n_flows);
    move_info.leaving_flow = LeavingFlowType("leaving_flow", n_flows);

    if (is_two_phase_flow)
    {
      this->mt_model = MassTransfer::MassTransferModel(
          MassTransfer::MTRType::Flowmap, liquid_scalar, gas_scalar);
    }
  }

  void SimulationUnit::setVolumes(std::span<const double> volumesgas,
                                  std::span<const double> volumesliq) const
  {

    std::span<double const> vg;
    this->liquid_scalar->setVolumes(volumesliq, flow_liquid->inverse_volume);
    if (gas_scalar)
    {
      this->gas_scalar->setVolumes(volumesgas, flow_gas->inverse_volume);
      vg = gas_scalar->getVolumeData();
    }
    else
    {
      vg = volumesgas;
    }

    std::span<double const> vl = liquid_scalar->getVolumeData();

    this->mc_unit->domain.setVolumes(vg, vl);
  }

  void SimulationUnit::post_init_concentration_file(const ScalarInitializer& scalar_init)
  {

    if (!scalar_init.liquid_buffer.has_value())
    {
      throw SimulationException(ErrorCodes::BadInitialiser);
    }

    if (!this->liquid_scalar->deep_copy_concentration(*scalar_init.liquid_buffer))
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
    flow_liquid = nullptr;
    flow_gas = nullptr;
    move_info = KernelInline::MoveInfo<ComputeSpace>();
  }

  std::size_t SimulationUnit::counter() const
  {
    return internal_counter_dead;
  }

  void SimulationUnit::set_probes(ProbeAutogeneratedBuffer&& _probes)
  {
    probes = std::move(_probes);
  }

  bool SimulationUnit::two_phase_flow() const
  {
    return is_two_phase_flow;
  }

  void SimulationUnit::post_init_concentration_functor(const ScalarInitializer& scalar_init)
  {
    // CmaRead::L2DView<double> cliq = this->liquid_scalar->getConcentrationView();

    auto cliqdata = this->liquid_scalar->getConcentrationData();

    auto get_view = [this](auto&& cdata)
    {
      return CmaRead::L2DView<double>({cdata.data(),
                                     static_cast<size_t>(cdata.size())},
                                    this->liquid_scalar->n_row(),
                                    this->liquid_scalar->n_col(),
                                    false);
    };

    auto cliq = get_view(cliqdata);

    CmaRead::L2DView<double> c;
    if (is_two_phase_flow)
    {
      assert(this->gas_scalar != nullptr);
      c = get_view(this->gas_scalar->getConcentrationData());
      assert(c.size() != 0);
    }

    for (size_t i = 0; i < cliq.getNCol(); ++i)
    {
      scalar_init.liquid_f_init.value()(i, cliq);
      if (is_two_phase_flow)
      {
        scalar_init.gas_f_init.value()(i, c);
      }
    }
  }

  void SimulationUnit::post_init_concentration(const ScalarInitializer& scalar_init)
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
    auto _compute_concentration = liquid_scalar->get_device_concentration();
    auto _containers = mc_unit->domain.data();

    Kokkos::parallel_for(
        "post_init_compartments",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,
                                                           mc_unit->domain.getNumberCompartments()),
        KOKKOS_LAMBDA(const int i) {
          // auto s = Kokkos::subview(_compute_concentration, i, Kokkos::ALL);
          _containers(i).concentrations = Kokkos::subview(_compute_concentration, Kokkos::ALL, i);
        });
    Kokkos::fence();
  }

  SimulationUnit::~SimulationUnit() = default;

  DiagonalView<ComputeSpace> SimulationUnit::get_kernel_diagonal()
  {

    return flow_liquid->get_kernel_diagonal();
  }

  CumulativeProbabilityView<ComputeSpace> SimulationUnit::get_kernel_cumulative_proba()
  {
    auto& matrix = flow_liquid->cumulative_probability;

    CumulativeProbabilityView<HostSpace> rd(matrix.data(), matrix.rows(), matrix.cols());

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), rd);
  }

  kernelContribution SimulationUnit::get_kernel_contribution()
  {
    return this->liquid_scalar->get_kernel_contribution();
  }

  void SimulationUnit::set_kernel_contribs_to_host()
  {
    this->liquid_scalar->set_kernel_contribs_to_host();
  }

  // void SimulationUnit::pimpl_deleter::operator()(ScalarSimulation* ptr) const
  // {
  //   delete ptr; // NOLINT
  // }

  [[nodiscard]] double& SimulationUnit::get_start_time_mut()
  {
    return starting_time;
  }
  [[nodiscard]] double& SimulationUnit::get_end_time_mut()
  {
    return end_time;
  }

  const Simulation::Feed::SimulationFeed& SimulationUnit::get_feed() const
  {
    return this->feed;
  }

} // namespace Simulation
