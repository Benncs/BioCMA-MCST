#include "common/thread_safe_container.hpp"
#include "mc/particles/particles_container.hpp"
#include "mc/unit.hpp"
#include "models/types.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <mc/domain.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng.hpp>
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

// TODO REMOVE
#include <iostream>
#include <utility>

namespace Simulation
{
  void pimpl_deleter::operator()(ScalarSimulation *ptr) const
  {
    delete ptr;
  }

  void initF(pimp_ptr_t &liq, pimp_ptr_t &gas)
  {
    bool host = gas != nullptr;
    for (auto i = 0; i < liq->C.cols(); ++i)
    {
      liq->C.coeffRef(1, i) = 8e-3;
      liq->C.coeffRef(0, static_cast<int>(i)) = 0.5; //0.5 G/L Glucose 
      if (host)
      {
        gas->C.coeffRef(1, i) = 321e-3 * 1e5 / (8.314 * (273.15 + 30)) * 0.2;
      }
    }

    liq->Mtot = liq->C * liq->getVolume();
    if (host)
    {
      gas->Mtot = gas->C * gas->getVolume();
    }
  }

  void p_kernel(double d_t,
                MC::MonteCarloUnit &unit,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::TheadSafeData> _extras,
                const KModel &_kmodel,
                const Simulation::move_kernel &_move_kernel,
                MC::Particles &p);

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        mc_container(std::move(other.mc_container)), host(other.host),
        n_thread(other.n_thread), flow_liquid(other.flow_liquid),
        flow_gas(other.flow_gas), kmodel(std::move(other.kmodel))
  {
  }

  SimulationUnit::SimulationUnit(
      const ExecInfo &info,
      std::unique_ptr<MC::MonteCarloUnit> &&_unit,
      std::unique_ptr<MC::ParticlesContainer> &&_container,
      size_t n_species,
      KModel _km,
      bool _host)
      : mc_unit(std::move(_unit)), mc_container(std::move(_container)),
        host(_host), n_thread(info.thread_per_process), flow_liquid(nullptr),
        flow_gas(nullptr), kmodel(std::move(_km))
  {
    if (this->mc_container->extras.empty())
    {
      throw std::runtime_error("Extra particle not initialised");
    }

    this->liquid_scalar =
        std::unique_ptr<ScalarSimulation, pimpl_deleter>(makeScalarSimulation(
            mc_unit->domain.getNumberCompartments(), n_species, n_thread));

    this->gas_scalar =
        (host) ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(
                     makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                                          n_species,
                                          0)) // No contribs for gas
               : nullptr;

    // FIXME
    initF(liquid_scalar, gas_scalar);

    post_init_container();
    post_init_compartments();
  }

  void SimulationUnit::setVolumes(std::span<double> volumesgas,
                                  std::span<double> volumesliq)
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

  void SimulationUnit::post_init_compartments()
  {
    int i = 0;
    size_t n_species = this->liquid_scalar->n_species();
    std::for_each(mc_unit->domain.begin(),
                  mc_unit->domain.end(),
                  [this, n_species, &i](auto &&cs)
                  {
                    cs.concentrations = std::span(
                        this->liquid_scalar->C.col(i).data(), n_species);
                    ++i;
                  });
  }

  void SimulationUnit::post_init_container()
  {
    auto &to_process = mc_container->to_process;
    const size_t n_compartments = mc_unit->domain.getNumberCompartments();
#pragma omp parallel for default(none) shared(to_process, n_compartments)      \
    num_threads(this->n_thread)
    for (auto it = to_process.begin(); it < to_process.end(); ++it)
    {
      auto &&particle = *it;
      particle.current_container =
          MC::uniform_int_rand(static_cast<size_t>(0), n_compartments - 1);

      auto &i_container = mc_unit->domain[particle.current_container];

      kmodel.init_kernel(particle);

      __ATOM_INCR__(i_container.n_cells);
    }
  }

  void SimulationUnit::cycleProcess(const double d_t)
  {

    const move_kernel _move_kernel = population_balance_flow(
        this->mc_unit->domain, this->flow_liquid);

    const auto &_kmodel = kmodel;

    auto &unit = *this->mc_unit;
    auto _contribs = this->liquid_scalar->getThreadContribs();
    auto &_extras = this->mc_container->extras;
    auto &to_process = mc_container->to_process;

#pragma omp parallel for num_threads(this->n_thread) default(none)             \
    shared(to_process, unit, _contribs, _extras, _kmodel, _move_kernel, d_t)   \
    schedule(static)
    for (auto it = to_process.begin(); it < to_process.end(); ++it)
    {
      p_kernel(d_t, unit, _contribs, _extras, _kmodel, _move_kernel, *it);
    }

    post_process_reducing();
  }

  void p_kernel(double d_t,
                MC::MonteCarloUnit &unit,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::TheadSafeData> _extras,
                const KModel &_kmodel,
                const move_kernel &_move_kernel,
                MC::Particles &p)
  {

    if (p.status == MC::CellStatus::DEAD)
    {
      return;
    }

    auto &domain = unit.domain;
    const double rnd = MC::double_unfiform();
    const double rdn2 = MC::double_unfiform();
    const size_t i_thread = omp_get_thread_num();
    auto &events = unit.ts_events[i_thread];
    auto &thread_contrib = _contribs[i_thread];
    auto &thread_extra = _extras[i_thread];
    const auto &concentrations = domain[p.current_container].concentrations;

    _move_kernel(rnd, rdn2, domain, p, d_t);

    _kmodel.update_kernel(d_t, p, concentrations);

    if (p.status == MC::CellStatus::DEAD)
    {
      events.incr<MC::EventType::Death>();
      __ATOM_DECR__(domain[p.current_container].n_cells)   // TODO: check overflow
      p.clearState(MC::CellStatus::DEAD);
      thread_extra.in_dead_state.emplace_back(&p);
    }
    else
    {
      if (p.status == MC::CellStatus::CYTOKINESIS)
      {
        events.incr<MC::EventType::NewParticle>();
        MC::Particles child = _kmodel.division_kernel(p);
        _kmodel.contribution_kernel(child, thread_contrib);
        __ATOM_INCR__(domain[child.current_container].n_cells)
        thread_extra.extra_process.emplace_back(std::move(child));
      }

      _kmodel.contribution_kernel(p, thread_contrib);
    }
  };

} // namespace Simulation
