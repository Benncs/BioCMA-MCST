#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <simulation/simulation.hpp>

#include <execution>
#include <iostream>
#include <omp.h>

#include <random>

SimulationUnit::SimulationUnit(
    std::unique_ptr<MC::MonteCarloUnit> &&_unit,
    std::unique_ptr<MC::ParticlesContainer> &&_container,
    bool _master)
    : unit(std::move(_unit)), container(std::move(_container)), gen(rd()),
      master(_master)
{

  // double XGlobal = 1.; // FIXME
  // double mtot = XGlobal * this->unit->domain.total_volume();
  // double m_particle = mtot / this->container->to_process.size();

  // for (auto &&particle : this->container->to_process)
  // {
  //   auto &container = unit->domain[particle.current_container];
  //   container.n_cells += 1;
  // }

  auto exec = std::execution::seq;
  std::for_each(exec,
                container->to_process.begin(),
                container->to_process.end(),
                [this](auto &&particle)
                {
                  auto &i_container = unit->domain[particle.current_container];
                  i_container.n_cells += 1;
                });
}

void SimulationUnit::cycle_process()
{

  size_t i = 0;
  std::uniform_real_distribution<> dis(0, 1.0);

  const auto krnl = [&dis, _gen = &this->gen](auto &&p)
  {
    auto rnd = dis(*_gen);
    if (rnd < 0.3)
    {
      p.id = 0;
    }
    else
    {
      p.id = 1;
    }
  };

  // #pragma omp parallel for
  //   for (auto it = container->to_process.begin();
  //        it < container->to_process.end();
  //        ++it)
  //   {
  //     krnl(*it);
  //   }

  auto exec = std::execution::par_unseq;
  std::for_each(
      exec, container->to_process.begin(), container->to_process.end(), krnl);
}
