#include "mc/particles/mcparticles.hpp"
#include "simulation/scalar_simulation.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <mc/prng.hpp>
#include <memory>
#include <omp.h>
#include <simulation/simulation.hpp>

// double XGlobal = 1.; // FIXME
// double mtot = XGlobal * this->unit->domain.total_volume();
// double m_particle = mtot / this->container->to_process.size();

// for (auto &&particle : this->container->to_process)
// {
//   auto &container = unit->domain[particle.current_container];
//   container.n_cells += 1;
// }

// auto exec = std::execution::seq;
// std::for_each(exec,
//               container->to_process.begin(),
//               container->to_process.end(),
//               [this](auto &&particle)
//               {
//                 particle.current_container = MC::intRand(size_t(0),
//                 unit->domain.n_compartments()-1); auto &i_container =
//                 unit->domain[particle.current_container];
//                 i_container.n_cells += 1;
//               });

// static void default_wrapper(MC::Particles &d)
// {
// }

// static void default_wrapper_p(MC::Particles &d, std::span<double const> dd)
// {
// }
namespace Simulation
{

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
  {
    unit = std::move(other.unit);
    container = std::move(other.container);
    host = other.host;
    flow_liquid = std::move(other.flow_liquid);
    flow_gas = std::move(other.flow_gas);
    kmodel = other.kmodel;
  }

  SimulationUnit::SimulationUnit(
      size_t n_species,
      std::unique_ptr<MC::MonteCarloUnit> &&_unit,
      std::unique_ptr<MC::ParticlesContainer> &&_container,
      const ExecInfo &info,
      bool _host)
      : unit(std::move(_unit)), container(std::move(_container)), host(_host)
  {

    this->liquid_scalar = std::make_unique<ScalarSimulation>(
        unit->domain.n_compartments(), n_species);
    np = info.thread_per_process;

    for (auto i = 0; i < liquid_scalar->C.cols(); ++i)
    {
      liquid_scalar->C.coeffRef(1, i) = 8e-3;
    }

    this->contribs =
        std::vector(info.thread_per_process,
                    Eigen::MatrixXd(this->liquid_scalar->C.rows(),
                                    this->liquid_scalar->C.cols()));
    for (auto &&i : contribs)
    {
      i.setZero();
    }

    if (host)
    {
      this->gas_scalar = std::make_unique<ScalarSimulation>(
          unit->domain.n_compartments(), n_species);
    }
    else
    {
      this->gas_scalar = nullptr;
    }
  }

  void SimulationUnit::post_init(KModel &&_km)
  {
    kmodel = std::move(_km);
    post_init_container();
    post_init_compartments();
  }

  void SimulationUnit::post_init_compartments()
  {
    int i = 0;
    size_t n_species = this->liquid_scalar->n_species();
    std::for_each(unit->domain.begin(),
                  unit->domain.end(),
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
          MC::uniform_int_rand(size_t(0), unit->domain.n_compartments() - 1);

      kmodel.init_kernel(particle);
      auto &i_container = unit->domain[particle.current_container];

#pragma omp atomic
      i_container.n_cells += 1;
    }
  }

  void SimulationUnit::cycle_process(const double d_t)
  {

    const move_kernel _move_kernel =
        pbf(*this->unit, *this->container, flow_liquid);

    const auto &_kmodel = kmodel;

    const auto &domain = this->unit->domain;

    auto &_contribs = this->contribs;

    const auto krnl =
        [&_move_kernel, d_t, _kmodel, &domain, &_contribs](auto &&p)
    {
      const double rnd = MC::double_unfiform();
      const double rdn2 = MC::double_unfiform();
      const size_t i_contrib = omp_get_thread_num();

      auto &mat = _contribs.at(i_contrib);
      const auto &concentrations = domain[p.current_container].concentrations;

      _move_kernel(rnd, rdn2, p, d_t);
      _kmodel.update_kernel(d_t, p, concentrations);
      _kmodel.division_kernel(p);
      _kmodel.contribution_kernel(p, mat);
    };

#pragma omp parallel for
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      krnl(*it);
    }

    for (auto &&i : _contribs)
    {
      this->liquid_scalar->biomass_contribution += i;
    }
  }

  void SimulationUnit::step(double d_t)
  {
    
    
    this->liquid_scalar->performStep(d_t, flow_liquid.transition_matrix);
    // this->gas_scalar->performStep(d_t, flow_gas.transition_matrix);
  }
} // namespace Simulation