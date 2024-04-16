#include "mc/particles/particles_container.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <mc/domain.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng.hpp>
#include <memory>
#include <omp.h>
#include <scalar_simulation.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transport.hpp>

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
      liq->C.coeffRef(0, static_cast<int>(i)) = 0.01;
      if (host)
      {
        gas->C.coeffRef(1, i) = 321e-3 * 1e5 / (8.314 * (273.15 + 30)) * 0.2;
      }
    }

    liq->Mtot = liq->C * liq->V;
    if (host)
    {
      gas->Mtot = gas->C * gas->V;
    }
  }

  void p_kernel(MC::ReactorDomain &domain,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::TheadSafeData> _extras,
                double d_t,
                const KModel &_kmodel,
                const Simulation::move_kernel &_move_kernel,
                MC::Particles &p);

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        container(std::move(other.container)), host(other.host),
        flow_liquid(std::move(other.flow_liquid)),
        flow_gas(std::move(other.flow_gas)), kmodel(other.kmodel)
  {
  }

  SimulationUnit::SimulationUnit(
      const ExecInfo &info,
      std::unique_ptr<MC::MonteCarloUnit> &&_unit,
      std::unique_ptr<MC::ParticlesContainer> &&_container,
      size_t n_species,
      bool _host)
      : mc_unit(std::move(_unit)), container(std::move(_container)),
        host(_host), n_thread(info.thread_per_process)
  {

    this->liquid_scalar =
        std::unique_ptr<ScalarSimulation, pimpl_deleter>(makeScalarSimulation(
            mc_unit->domain.n_compartments(), n_species, n_thread));

    this->gas_scalar =
        (host) ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(
                     makeScalarSimulation(mc_unit->domain.n_compartments(),
                                          n_species,
                                          0)) // No contribs for gas
               : nullptr;

    // FIXME
    initF(liquid_scalar, gas_scalar);
  }

  void SimulationUnit::postInit(KModel &&_km)
  {
    kmodel = std::move(_km);
    post_init_container();
    post_init_compartments();
  }

  void SimulationUnit::setVolumes(std::vector<double> &&volumesgas,
                                  std::vector<double> &&volumesliq)
  {

    std::span<double> vg;
    this->liquid_scalar->setV(std::move(volumesliq));
    if (gas_scalar)
    {
      this->gas_scalar->setV(std::move(volumesgas));
      vg = std::span<double>(gas_scalar->V.diagonal().data(),
                             this->mc_unit->domain.n_compartments());
    }
    else
    {
      vg = volumesgas;
    }

    std::span<double> vl =
        std::span<double>(liquid_scalar->V.diagonal().data(),
                          this->mc_unit->domain.n_compartments());

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

#pragma omp parallel for
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      auto &&particle = *it;
      particle.current_container =
          MC::uniform_int_rand(size_t(0), mc_unit->domain.n_compartments() - 1);

      auto &i_container = mc_unit->domain[particle.current_container];

      kmodel.init_kernel(particle);

#pragma omp atomic
      i_container.n_cells += 1;
    }
  }

  void SimulationUnit::cycleProcess(const double d_t)
  {

    const move_kernel _move_kernel =
        pbf(*this->mc_unit, *this->container, flow_liquid);

    const auto &_kmodel = kmodel;

    auto &domain = this->mc_unit->domain;
    auto &_contribs = this->liquid_scalar->contribs;
    auto &_extras = this->container->extras;

    const auto krnl =
        [&_move_kernel, d_t, _kmodel, &domain, &_contribs, &_extras](
            auto &&p) -> void
    { p_kernel(domain, _contribs, _extras, d_t, _kmodel, _move_kernel, p); };

    execute_process_knrl(krnl);
    post_process_reducing();
  }

  void p_kernel(MC::ReactorDomain &domain,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::TheadSafeData> _extras,
                double d_t,
                const KModel &_kmodel,
                const move_kernel &_move_kernel,
                MC::Particles &p)
  {

    if (p.status == MC::CellStatus::DEAD)
    {
      return;
    }

    const double rnd = MC::double_unfiform();
    const double rdn2 = MC::double_unfiform();
    const size_t i_thread = omp_get_thread_num();

    auto &thread_contrib = _contribs[i_thread];
    auto &thread_extra = _extras[i_thread];
    const auto &concentrations = domain[p.current_container].concentrations;

    _move_kernel(rnd, rdn2, p, d_t);

    _kmodel.update_kernel(d_t, p, concentrations);

    if (p.status != MC::CellStatus::DEAD)
    {
      if (p.status == MC::CellStatus::CYTOKINESIS)
      {

        MC::Particles child = _kmodel.division_kernel(p);

        _kmodel.contribution_kernel(child, thread_contrib);

#pragma omp atomic
        domain[child.current_container].n_cells += 1;

        thread_extra.extra_process.emplace_back(std::move(child));
      }

      _kmodel.contribution_kernel(p, thread_contrib);
    }
    else
    {
      p.clearState(MC::CellStatus::DEAD);
      thread_extra.in_dead_state.emplace_back(&p);

// TODO CHECK POSSIBLE OVERFLOW
#pragma omp atomic
      domain[p.current_container].n_cells -= 1;
    }
  };

} // namespace Simulation