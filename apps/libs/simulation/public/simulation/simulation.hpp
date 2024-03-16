#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <cstdlib>
#include <memory>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>

#include <random>

class SimulationUnit
{
public:
  explicit SimulationUnit(std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                          std::unique_ptr<MC::ParticlesContainer> &&_container,
                          bool master = false);
  std::unique_ptr<MC::MonteCarloUnit> unit;
  std::unique_ptr<MC::ParticlesContainer> container;
  SimulationUnit(SimulationUnit &&other) noexcept
  {
    unit = std::move(other.unit);
    container = std::move(other.container);
    gen = other.gen;
    master = other.master;
  }
  void cycle_process();
  int rng = rand() % ((10 + 1) - 1) + 1;
  std::random_device rd;
  std::mt19937 gen;

private:
  bool master;
};

#endif //__SIMULATIONS_UNIT_HPP__