#include "cma_read/light_2d_view.hpp"
#include "common/thread_safe_container.hpp"
#include "mc/prng/prng.hpp"
#include "models/types.hpp"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execution>
#include <mc/domain.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/thread_private_data.hpp>
#include <mc/unit.hpp>

#include <memory>

#include <scalar_simulation.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>
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

      liq(0, i) = 9;
    }
    else
    {
      liq(0, i) = 5.;
    }
  }

  void init_f_mat_gas(size_t i, CmaRead::L2DView<double> &gas)
  {
    gas(1, i) = 15e-3 * 1e5 / (8.314 * (273.15 + 30)) * 0.2;
  }

  void p_kernel(size_t current_i_particle,
                double d_t,
                MC::ReactorDomain &domain,
                std::span<MC::EventContainer> ts_events,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::ThreadPrivateData> _extras,
                const KModel &_kmodels,
                MC::Particles &p,
                   const auto& diag,
                auto &cumulative_probability);

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        is_two_phase_flow(other.is_two_phase_flow), n_thread(other.n_thread),
        flow_liquid(other.flow_liquid), flow_gas(other.flow_gas),
        kmodel(std::move(other.kmodel))
  {
  }

  SimulationUnit::SimulationUnit(const ExecInfo &info,
                                 std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                                 std::span<double> volumesgas,
                                 std::span<double> volumesliq,
                                 size_t n_species,
                                 KModel _km,
                                 bool _gas_flow)
      : mc_unit(std::move(_unit)), is_two_phase_flow(_gas_flow),
        n_thread(info.thread_per_process), flow_liquid(nullptr),
        flow_gas(nullptr), kmodel(std::move(_km))
  {

    if (this->mc_unit->extras.empty())
    {
      throw std::runtime_error("Extra particle not initialised");
    }

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

  void SimulationUnit::post_init_concentration()
  {
    auto rng = MC::PRNG();
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
    int i = 0;
    size_t n_species = this->liquid_scalar->n_species();
    std::for_each(mc_unit->domain.begin(),
                  mc_unit->domain.end(),
                  [this, n_species, &i](auto &&cs)
                  {
                    cs.concentrations = std::span(
                        this->liquid_scalar->concentration.col(i).data(),
                        n_species);
                    ++i;
                  });
  }

  void SimulationUnit::cycleProcess(const double d_t)
  {
    auto &to_process = this->mc_unit->container.to_process;

 

    auto view_cumulative_probability = this->flow_liquid->get_view_cum_prob();


    

    const auto& diag = this->flow_liquid->diag_transition;
    auto contribs = this->liquid_scalar->getThreadContribs();
    auto &extras = this->mc_unit->extras;
    auto &ts_events = this->mc_unit->ts_events;
    const auto size = static_cast<size_t>(to_process.size());
    auto &domain = this->mc_unit->domain;

#pragma omp for
    for (size_t i_particle = 0; i_particle < size; ++i_particle)
    {
      p_kernel(i_particle,
               d_t,
               domain,
               ts_events,
               contribs,
               extras,
               kmodel,
               to_process[i_particle],
               diag,
               view_cumulative_probability);
    }


#pragma omp master
    {
      post_process_reducing();
    }
  }

  void p_kernel(size_t current_i_particle,
                double d_t,
                MC::ReactorDomain &domain,
                std::span<MC::EventContainer> ts_events,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::ThreadPrivateData> _extras,
                const KModel &_kmodel,
                MC::Particles &particle,
                const auto& diag,
                auto &cumulative_probability)
  {

    if (particle.status == MC::CellStatus::DEAD)
    {
      return;
    }

    const size_t i_thread = omp_get_thread_num();
    auto &events = ts_events[i_thread];
    auto &thread_contrib = _contribs[i_thread];
    auto &thread_extra = _extras[i_thread];

    auto &rng = thread_extra.rng;
    const double random_number_1 = rng.double_unfiform();
    const double random_number_2 = rng.double_unfiform();
    const double random_number_3 = rng.double_unfiform();

    auto register_spawned_particle = [&](MC::Particles &&p)
    {
      events.incr<MC::EventType::NewParticle>();
      _kmodel.contribution_kernel(p, thread_contrib);
      __ATOM_INCR__(domain[p.current_container].n_cells)
      thread_extra.extra_process.emplace_back(std::move(p));
    };

    kernel_move(random_number_1,
                random_number_2,
                domain,
                particle,
                d_t,
                diag,
                cumulative_probability);

    kernel_exit(d_t, random_number_3, domain, particle);

    if (particle.status == MC::CellStatus::OUT)
    {
      events.incr<MC::EventType::Exit>();
      __ATOM_DECR__(
          domain[particle.current_container].n_cells) // TODO: check overflow
      particle.clearState(MC::CellStatus::DEAD);
      thread_extra.index_in_dead_state.emplace_back(current_i_particle);
      return;
    }

    _kmodel.update_kernel(
        d_t, particle, domain[particle.current_container].concentrations);

    if (particle.status == MC::CellStatus::DEAD)
    {
      events.incr<MC::EventType::Death>();
      __ATOM_DECR__(
          domain[particle.current_container].n_cells) // TODO: check overflow

      particle.clearState(MC::CellStatus::DEAD);
      thread_extra.index_in_dead_state.emplace_back(current_i_particle);
    }
    else
    {
      if (particle.status == MC::CellStatus::CYTOKINESIS)
      {
        particle.status = MC::CellStatus::IDLE;
        register_spawned_particle(_kmodel.division_kernel(particle));
      }

      _kmodel.contribution_kernel(particle, thread_contrib);
    }
  };

  void SimulationUnit::pimpl_deleter::operator()(ScalarSimulation *ptr) const
  {
    delete ptr;
  }

} // namespace Simulation
