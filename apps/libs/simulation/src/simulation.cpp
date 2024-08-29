#include "common/kokkos_vector.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Random.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>
#include <mc/thread_private_data.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <scalar_simulation.hpp>
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

  void init_f_mat_liq(size_t i, CmaRead::L2DView<double> &liq)
  {
    if (i == 0)
    {

      liq(0, i) = 20;
    }
    else
    {
      liq(0, i) = 0.1;
    }
  }

  void init_f_mat_gas(size_t i, CmaRead::L2DView<double> &gas)
  {
    gas(1, i) = 15e-3 * 1e5 / (8.314 * (273.15 + 30)) * 0.2;
  }

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        is_two_phase_flow(other.is_two_phase_flow), n_thread(other.n_thread),
        flow_liquid(other.flow_liquid), flow_gas(other.flow_gas)
  {
  }

  SimulationUnit::SimulationUnit(const ExecInfo &info,
                                 std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                                 std::span<double> volumesgas,
                                 std::span<double> volumesliq,
                                 size_t n_species,
                                 bool _gas_flow)
      : mc_unit(std::move(_unit)), is_two_phase_flow(_gas_flow),
        n_thread(info.thread_per_process), flow_liquid(nullptr),
        flow_gas(nullptr)
  {

    this->liquid_scalar = std::unique_ptr<ScalarSimulation, pimpl_deleter>(
        makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                             n_species,
                             n_thread,
                             volumesliq));

    this->gas_scalar =
        (is_two_phase_flow)
            ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(
                  makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                                       n_species,
                                       0,
                                       volumesgas)) // No contribs for gas
            : nullptr;

    post_init_concentration();
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

  Kokkos::View<const size_t **, Kokkos::LayoutStride, ComputeSpace>
  SimulationUnit::get_view_neighbor() const
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

  void SimulationUnit::post_init_concentration()
  {
    // auto rng = MC::PRNG();
    CmaRead::L2DView<double> cliq = this->liquid_scalar->getConcentrationView();
    CmaRead::L2DView<double> *cgas = nullptr;
    if (is_two_phase_flow)
    {
      *cgas = this->gas_scalar->getConcentrationView();
    }

    for (size_t i = 0; i < cliq.getNCol(); ++i)
    {
      init_f_mat_liq(i, cliq);
      if (is_two_phase_flow)
      {
        init_f_mat_gas(i, *cgas);
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

  Kokkos::View<double *,
               Kokkos::LayoutLeft,
               HostSpace,
               Kokkos::MemoryTraits<Kokkos::RandomAccess>>
  SimulationUnit::get_dg()
  {
    Kokkos::View<double *,
                 Kokkos::LayoutLeft,
                 HostSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        diag(this->flow_liquid->diag_transition.data(),
             this->flow_liquid->diag_transition.size());

    diag_transition = Kokkos::create_mirror_view_and_copy(ComputeSpace(), diag);

    return diag_transition;
  }

  Kokkos::View<double **,
               Kokkos::LayoutLeft,
               HostSpace,
               Kokkos::MemoryTraits<Kokkos::RandomAccess>>
  SimulationUnit::get_cp()
  {
    auto &matrix = flow_liquid->cumulative_probability;

    Kokkos::View<double **,
                 Kokkos::LayoutLeft,
                 HostSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        rd(matrix.data(), matrix.rows(), matrix.cols());

    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), rd);
  }

  Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace>
  SimulationUnit::get_contribs()
  {
    return Kokkos::create_mirror_view_and_copy(ComputeSpace(),
                                               liquid_scalar->k_contribs);
  }

  void SimulationUnit::set_contribs(
      Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace> c)
  {
    Kokkos::deep_copy(liquid_scalar->k_contribs, c);
  }

  void SimulationUnit::pimpl_deleter::operator()(ScalarSimulation *ptr) const
  {
    delete ptr;
  }

} // namespace Simulation
