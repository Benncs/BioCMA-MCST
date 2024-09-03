#include "simulation/simulation_exception.hpp"
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
#include <scalar_simulation.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>
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

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        is_two_phase_flow(other.is_two_phase_flow),
        flow_liquid(other.flow_liquid), flow_gas(other.flow_gas)
  {
  }

  SimulationUnit::SimulationUnit(const ExecInfo &info,
                                 std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                                 ScalarInitializer scalar_init)
      : mc_unit(std::move(_unit)), is_two_phase_flow(scalar_init.gas_flow),

        flow_liquid(nullptr), flow_gas(nullptr)
  {
    this->liquid_scalar = std::unique_ptr<ScalarSimulation, pimpl_deleter>(
        makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                             scalar_init.n_species,
                             scalar_init.volumesliq));

    this->gas_scalar =
        (is_two_phase_flow)
            ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(
                  makeScalarSimulation(
                      mc_unit->domain.getNumberCompartments(),
                      scalar_init.n_species,
                      scalar_init.volumesgas)) // No contribs for gas
            : nullptr;

    post_init_concentration(scalar_init.liquid_f_init, scalar_init.gaz_f_init);
    post_init_compartments();
    const std::size_t n_exit_flow = 1;
    index_leaving_flow =
        Kokkos::View<size_t *, ComputeSpace>("index_leaving_flow", n_exit_flow);

    leaving_flow =
        Kokkos::View<double *, ComputeSpace>("leaving_flow", n_exit_flow);
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

  SimulationUnit::NeighborsViewCompute
  SimulationUnit::get_kernel_neighbors() const
  {
    const auto view_neighbors = this->mc_unit->domain.getNeighbors();

    const Kokkos::LayoutStride layout(view_neighbors.getNRow(),
                                      view_neighbors.getNCol(),
                                      view_neighbors.getNCol(),
                                      1);

    const auto host_view = Kokkos::View<const size_t **, Kokkos::LayoutStride>(
        view_neighbors.data().data(), layout);

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_view);
  }

  void SimulationUnit::post_init_concentration(init_scalar_f_t liquid,
                                               init_scalar_f_t gas)
  {
    CmaRead::L2DView<double> cliq = this->liquid_scalar->getConcentrationView();
    CmaRead::L2DView<double> *cgas = nullptr;
    if (is_two_phase_flow)
    {
      *cgas = this->gas_scalar->getConcentrationView();
    }

    for (size_t i = 0; i < cliq.getNCol(); ++i)
    {
      liquid(i, cliq);
      if (is_two_phase_flow)
      {
        gas(i, *cgas);
      }
    }

    if ((this->liquid_scalar->concentration.array() < 0).any())
    {
      throw SimulationException(ErrorCodes::BadConcentrationInitLiq);
    }

    if (is_two_phase_flow)
    {
      if ((this->gas_scalar->concentration.array() < 0).any())
      {
        throw SimulationException(ErrorCodes::BadConcentrationInitGas);
      }
    }

    this->liquid_scalar->total_mass =
        this->liquid_scalar->concentration * this->liquid_scalar->getVolume();
    if (is_two_phase_flow)
    {
      this->gas_scalar->total_mass =
          this->gas_scalar->concentration * this->gas_scalar->getVolume();
    }
  }

  void SimulationUnit::post_init_compartments()
  {
    auto _compute_concentration = liquid_scalar->compute_concentration;
    auto _containers = mc_unit->domain.data();

    Kokkos::parallel_for(
        "post_init_compartments",
        mc_unit->domain.getNumberCompartments(),
        KOKKOS_LAMBDA(const int i) {
          _containers(i).concentrations =
              Kokkos::subview(_compute_concentration, i, Kokkos::ALL);
        });
    Kokkos::fence();
    // int i = 0;
    // size_t n_species = this->liquid_scalar->n_species();
    // std::for_each(mc_unit->domain.begin(),
    //               mc_unit->domain.end(),
    //               [this, n_species, &i](auto &&cs)
    //               {
    //                 cs.concentrations = std::span(
    //                     this->liquid_scalar->concentration.col(i).data(),
    //                     n_species);
    //                 ++i;
    //               });
  }

  SimulationUnit::DiagonalViewCompute SimulationUnit::get_kernel_diagonal()
  {
    DiagonalViewCompute diag(this->flow_liquid->diag_transition.data(),
                             this->flow_liquid->diag_transition.size());

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), diag);
  }

  SimulationUnit::CumulativeProbabilityViewCompute
  SimulationUnit::get_kernel_cumulative_proba()
  {
    auto &matrix = flow_liquid->cumulative_probability;

    CumulativeProbabilityViewCompute rd(
        matrix.data(), matrix.rows(), matrix.cols());

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), rd);
  }

  SimulationUnit::ContributionViewCompute
  SimulationUnit::get_kernel_contribution()
  {
    return Kokkos::create_mirror_view_and_copy(ComputeSpace(),
                                               liquid_scalar->k_contribs);
  }

  void SimulationUnit::set_kernel_contribs_to_host(ContributionViewCompute c)
  {
    Kokkos::deep_copy(liquid_scalar->k_contribs, c);
  }

  void SimulationUnit::pimpl_deleter::operator()(ScalarSimulation *ptr) const
  {
    delete ptr;
  }

} // namespace Simulation
