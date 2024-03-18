#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <random>
#include <simulation/transport.hpp>

#include <simulation/scalar_simulation.hpp>

namespace Simulation
{
  class SimulationUnit
  {
  public:
    explicit SimulationUnit(size_t n_species,
        std::unique_ptr<MC::MonteCarloUnit> &&_unit,
        std::unique_ptr<MC::ParticlesContainer> &&_container,
        bool host = false);

    std::unique_ptr<MC::MonteCarloUnit> unit;
    std::unique_ptr<MC::ParticlesContainer> container;
    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;

    auto getC(){return this->_implScalar->C;}
   
    void step(double d_t);
    void cycle_process(double d_t);
    std::random_device rd;
    std::mt19937 gen;

    void setFLows(MatFlow *_flows)
    {
      flows = _flows;
    }

  private:
    bool host;
    MatFlow *flows; // TODO OPTI
    std::unique_ptr<ScalarSimulation> _implScalar;
  };
} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__