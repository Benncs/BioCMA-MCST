#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <models/types.hpp>

namespace Simulation
{
  class MatFlow;
  class ScalarSimulation;
  struct pimpl_deleter
  {
    void operator()(ScalarSimulation *) const;
  };

  using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;

  void initF(pimp_ptr_t &liq, pimp_ptr_t &gas);

  class SimulationUnit
  {
  public:
    explicit SimulationUnit(
        const ExecInfo &info,
        std::unique_ptr<MC::MonteCarloUnit> &&_unit,
        std::unique_ptr<MC::ParticlesContainer> &&_container,
        size_t n_species,
        bool host = false);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;
    std::unique_ptr<MC::ParticlesContainer> container;
    ReactorState *state = nullptr;

    void postInit(KModel &&_km);

    std::span<double> getCliqData();

    void setVolumes(std::vector<double> &&volumesgas,
                    std::vector<double> &&volumesliq);

    void step(double d_t);

    void cycleProcess(double d_t);

    void setLiquidFlow(MatFlow *_flows_l);

    void setGasFlow(MatFlow *_flows_g);

    std::span<double> get_contributionData();

    void reduceContribs(std::span<double> data, size_t n_rank);

    void clearContribution();

  private:
    void post_init_container();
    void post_init_compartments();

    void post_process_reducing();
    void execute_process_knrl(const auto &kernel);

    bool host;
    size_t n_thread;

    MatFlow *flow_liquid; // TODO OPTI
    MatFlow *flow_gas;    // TODO OPTI
    KModel kmodel;

    std::unique_ptr<ScalarSimulation, pimpl_deleter> liquid_scalar;
    std::unique_ptr<ScalarSimulation, pimpl_deleter> gas_scalar;
  };

  inline void SimulationUnit::setLiquidFlow(MatFlow *_flows_l)
  {
    flow_liquid = _flows_l;
  }

  inline void SimulationUnit::setGasFlow(MatFlow *_flows_g)
  {
    flow_gas = _flows_g;
  }

  inline void SimulationUnit::execute_process_knrl(const auto &kernel)
  {
#pragma omp parallel for num_threads(n_thread)
    for (auto it = container->to_process.begin();
         it < container->to_process.end();
         ++it)
    {
      kernel(*it);
    }
  }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__