#include "simulation/feed_descriptor.hpp"
#include "simulation/scalar_initializer.hpp"
#include "simulation/simulation_exception.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Printf.hpp>
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
#include <simulation/simulation.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
#include <transport.hpp>

#ifndef USE_PYTHON_MODULE
#  include <omp.h>
#else
#  define omp_get_thread_num() 0 // 1 thread
#endif

#include <utility>

namespace Simulation
{

  SimulationUnit::SimulationUnit(SimulationUnit&& other) noexcept
      : mc_unit(std::move(other.mc_unit)), flow_liquid(other.flow_liquid), flow_gas(other.flow_gas),
        is_two_phase_flow(other.is_two_phase_flow), internal_counter_dead(std::move(other.internal_counter_dead)),
        waiting_allocation_particle(std::move(other.waiting_allocation_particle))
  {
  }

  SimulationUnit::SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                                 const ScalarInitializer& scalar_init,
                                 std::optional<Feed::SimulationFeed> _feed)
      : mc_unit(std::move(_unit)), flow_liquid(nullptr),

        flow_gas(nullptr), is_two_phase_flow(scalar_init.gas_flow),internal_counter_dead("internal_counter_dead"),
        waiting_allocation_particle("waiting_allocation_particle")
  {
    this->liquid_scalar = std::unique_ptr<ScalarSimulation, pimpl_deleter>(makeScalarSimulation(
        mc_unit->domain.getNumberCompartments(), scalar_init.n_species, scalar_init.volumesliq));

    this->gas_scalar = (is_two_phase_flow)
                           ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(makeScalarSimulation(
                                 mc_unit->domain.getNumberCompartments(),
                                 scalar_init.n_species,
                                 scalar_init.volumesgas)) // No contribs for gas
                           : nullptr;

    post_init_concentration(scalar_init);
    post_init_compartments();

    this->feed =
        _feed.has_value() ? std::move(*_feed) : Feed::SimulationFeed{std::nullopt, std::nullopt};

    size_t n_flows = 0;
    if (this->feed.liquid.has_value())
    {
      n_flows = this->feed.liquid->size();
    }
    index_leaving_flow = LeavingFlowIndexType("index_leaving_flow", n_flows);

    leaving_flow = LeavingFlowType("leaving_flow", n_flows);
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

  NeighborsViewCompute SimulationUnit::get_kernel_neighbors() const
  {
    const auto view_neighbors = this->mc_unit->domain.getNeighbors();

    const Kokkos::LayoutStride layout(
        view_neighbors.getNRow(), view_neighbors.getNCol(), view_neighbors.getNCol(), 1);

    // const Kokkos::LayoutLeft layout(view_neighbors.getNRow(),
    //                                   view_neighbors.getNCol());

    const auto host_view = NeighborsView<HostSpace>(view_neighbors.data().data(), layout);

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_view);
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

  void SimulationUnit::post_init_concentration_functor(const ScalarInitializer& scalar_init)
  {
    CmaRead::L2DView<double> cliq = this->liquid_scalar->getConcentrationView();
    CmaRead::L2DView<double> c;
    if (is_two_phase_flow)
    {
      assert(this->gas_scalar != nullptr);
      c = this->gas_scalar->getConcentrationView();
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

  DiagonalViewCompute SimulationUnit::get_kernel_diagonal()
  {

    return flow_liquid->get_kernel_diagonal();
  }

  CumulativeProbabilityViewCompute SimulationUnit::get_kernel_cumulative_proba()
  {
    auto& matrix = flow_liquid->cumulative_probability;

    CumulativeProbabilityView<HostSpace> rd(matrix.data(), matrix.rows(), matrix.cols());

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), rd);
  }

  kernelContribution SimulationUnit::get_kernel_contribution()
  {
    return this->liquid_scalar->get_kernel_contribution();
  }

  void SimulationUnit::set_kernel_contribs_to_host(kernelContribution c)
  {
    this->liquid_scalar->set_kernel_contribs_to_host(std::move(c));
  }

  void SimulationUnit::pimpl_deleter::operator()(ScalarSimulation* ptr) const
  {
    delete ptr; // NOLINT
  }

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
