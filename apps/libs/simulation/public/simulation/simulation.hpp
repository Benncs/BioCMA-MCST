#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include "mc/particles/mcparticles.hpp"
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <simulation/transport.hpp>
#include <scalar_simulation.hpp>
#include <simulation/models/types.hpp>

namespace Simulation
{

  class ScalarSimulation; 

  class SimulationUnit
  {
  public:
    explicit SimulationUnit(
        size_t n_species,
        std::unique_ptr<MC::MonteCarloUnit> &&_unit,
        std::unique_ptr<MC::ParticlesContainer> &&_container,
        const ExecInfo &info,
        bool host = false);

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;
    std::unique_ptr<MC::ParticlesContainer> container;

    void post_init(KModel &&_km);

    auto &getCliq();

    std::span<double> getCliqData();


    auto &getCgas();

    void setVolumes(std::vector<double> &&volumesgas,
                    std::vector<double> &&volumesliq);

    void step(double d_t);
    void cycle_process(const double d_t);

    void setLiquidFlow(MatFlow &&_flows_l);

    void setGasFlow(MatFlow &&_flows_g);

    auto &get_contribution();

    std::span<double>  get_contributionData();

    void reduce_contribs(std::span<double> data ,size_t n_rank);

    void clear_contribution();
    ReactorState *state = nullptr;

  private:
    void post_init_container();
    void post_init_compartments();

    void post_process_reducing();
    void execute_process_knrl(const auto &f);

    bool host;
    KModel kmodel;
    size_t np;

    MatFlow flow_liquid; // TODO OPTI
    MatFlow flow_gas;    // TODO OPTI
    std::unique_ptr<ScalarSimulation> liquid_scalar;
    std::unique_ptr<ScalarSimulation> gas_scalar;

    std::vector<Eigen::MatrixXd> contribs;
    std::vector<std::vector<MC::Particles>> extras_p;
  };

  inline void SimulationUnit::setLiquidFlow(MatFlow &&_flows_l)
  {
    flow_liquid = std::move(_flows_l);
  }

  inline void SimulationUnit::setGasFlow(MatFlow &&_flows_g)
  {
    flow_gas = std::move(_flows_g);
  }

  

  inline void SimulationUnit::execute_process_knrl(const auto &f)
  {
#pragma omp parallel for num_threads(this->contribs.size())
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      f(*it);
    }
  }



} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__