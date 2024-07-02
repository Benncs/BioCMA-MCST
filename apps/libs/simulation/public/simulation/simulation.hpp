#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include "mc/prng/prng.hpp"

#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <models/types.hpp>

namespace Simulation
{
  class PreCalculatedHydroState;
  class ScalarSimulation;
  struct pimpl_deleter
  {
    void operator()(ScalarSimulation *) const;
  };

  using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;

  void initF(pimp_ptr_t &liq, pimp_ptr_t &gas); // TODO Delete

  class SimulationUnit
  {
  public:
    explicit SimulationUnit(const ExecInfo &info,
                            std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                            std::span<double> volumesgas,
                            std::span<double> volumesliq,
                            size_t n_species,
                            KModel _km,MC::DistributionVariantInt&& param,
                            bool _gas_flow = false);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;
    SimulationUnit &operator=(SimulationUnit &&rhs) = delete;
    SimulationUnit &operator=(const SimulationUnit &rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] const KModel &getModel() const;

    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::tuple<size_t, size_t> getDim() const;
    [[nodiscard]] std::span<double> getCgasData() const;
    [[nodiscard]] std::span<double> getContributionData() const;

    void setVolumes(std::span<const double> volumesgas, std::span<const double> volumesliq);

    void step(double d_t, const CmaRead::ReactorState &state);

    void cycleProcess(double d_t);

    void setLiquidFlow(PreCalculatedHydroState *_flows_l);

    void setGasFlow(PreCalculatedHydroState *_flows_g);

    void reduceContribs(std::span<double> data, size_t n_rank);

    void clearContribution();

  private:
    void post_init_container(MC::DistributionVariantInt distribution_variant);
    void post_init_compartments();

    void post_process_reducing();
    void execute_process_knrl(const auto &kernel);

    bool is_two_phase_flow;
    size_t n_thread;

    PreCalculatedHydroState *flow_liquid; // TODO OPTI
    PreCalculatedHydroState *flow_gas;    // TODO OPTI
    KModel kmodel;

    std::vector<KModel> private_models;

    std::unique_ptr<ScalarSimulation, pimpl_deleter> liquid_scalar;
    std::unique_ptr<ScalarSimulation, pimpl_deleter> gas_scalar;
  };

  inline const KModel &SimulationUnit::getModel() const
  {
    return kmodel;
  }

  inline void SimulationUnit::setLiquidFlow(PreCalculatedHydroState *_flows_l)
  {
    flow_liquid = _flows_l;
  }

  inline void SimulationUnit::setGasFlow(PreCalculatedHydroState *_flows_g)
  {
    flow_gas = _flows_g;
  }

  inline void SimulationUnit::execute_process_knrl(const auto &kernel)
  {
    auto &container = mc_unit->container;
#pragma omp parallel for num_threads(n_thread) default(none)                   \
    shared(kernel, container) schedule(static)
    for (auto it = container.to_process.begin();
         it < container.to_process.end();
         ++it)
    {
      kernel(*it);
    }
  }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
