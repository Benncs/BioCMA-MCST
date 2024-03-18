#include "simulation/scalar_simulation.hpp"
#include <algorithm>
#include <memory>
#include <simulation/simulation.hpp>

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <random>

#include <iostream>
#include <mc/prng.hpp>

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

namespace Simulation
{

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
  {
    unit = std::move(other.unit);
    container = std::move(other.container);
    gen = other.gen;
    host = other.host;
  }

  SimulationUnit::SimulationUnit(
      size_t n_species,
      std::unique_ptr<MC::MonteCarloUnit> &&_unit,
      std::unique_ptr<MC::ParticlesContainer> &&_container,
      bool _host)
      : unit(std::move(_unit)), container(std::move(_container)), gen(rd()),
        host(_host)
  {
    this->_implScalar = std::make_unique<ScalarSimulation>(
        unit->domain.n_compartments(), n_species);

#pragma omp parallel for
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      auto &&particle = *it;
      particle.current_container =
          MC::intRand(size_t(0), unit->domain.n_compartments() - 1);
      auto &i_container = unit->domain[particle.current_container];
#pragma omp critical
      {
        i_container.n_cells += 1;
      }
    }

    int i = 0;

    std::for_each(unit->domain.begin(),
                  unit->domain.end(),
                  [this, n_species, &i](auto &&cs)
                  {
                    cs.concentrations = std::span(
                        this->_implScalar->C.col(i).data(), n_species);
                    ++i;
                  });
  }

  void SimulationUnit::cycle_process(double d_t)
  {

    size_t i = 0;

    move_kernel _move_kernel = pbf(*this->unit, *this->container, *flows);

    const auto krnl = [_gen = &this->gen, &_move_kernel, d_t](auto &&p)
    {
      std::uniform_real_distribution<> dis(0, 1.0);

      auto rnd = MC::doubleRand(0, 1);
      auto rdn2 = MC::doubleRand(0, 1);
      _move_kernel(rnd, rdn2, p, d_t);
    };

    // auto exec = std::execution::par_unseq;
    // std::for_each(
    //     exec, container->to_process.begin(), container->to_process.end(),
    //     krnl);

#pragma omp parallel for
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      krnl(*it);
    }
  }

  void SimulationUnit::step(double d_t)
  {
    Eigen::MatrixXd phi_v(this->_implScalar->C.rows(),
                          this->_implScalar->C.cols());

    phi_v.setZero();

    this->_implScalar->performStep(d_t, flows->transition_matrix, phi_v);
  }
} // namespace Simulation