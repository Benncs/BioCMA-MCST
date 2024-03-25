#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__



#include <reactorstate.hpp>
#include <common/common.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <simulation/transport.hpp>

#include <simulation/models/types.hpp>
#include <simulation/scalar_simulation.hpp>

namespace Simulation
{
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

    std::unique_ptr<MC::MonteCarloUnit> unit;
    std::unique_ptr<MC::ParticlesContainer> container;

    void post_init(KModel &&_km);

    auto &getCliq();

    auto &getCgas();

    void setVolumes(std::vector<double> &&volumesgas,
                    std::vector<double> &&volumesliq);

    void step(double d_t);
    void cycle_process(const double d_t);

    void setLiquidFlow(MatFlow &&_flows_l);

    void setGasFlow(MatFlow &&_flows_g);

    auto &get_contribution();

    void clear_contribution();
    ReactorState *state = nullptr;

  private:
    void post_init_container();
    void post_init_compartments();
    bool host;
    KModel kmodel;
    size_t np;

    MatFlow flow_liquid; // TODO OPTI
    MatFlow flow_gas;    // TODO OPTI
    std::unique_ptr<ScalarSimulation> liquid_scalar;
    std::unique_ptr<ScalarSimulation> gas_scalar;

    std::vector<Eigen::MatrixXd> contribs;
  };

  inline void SimulationUnit::setLiquidFlow(MatFlow &&_flows_l)
  {
    flow_liquid = std::move(_flows_l);
  }

  inline void SimulationUnit::setGasFlow(MatFlow &&_flows_g)
  {
    flow_gas = std::move(_flows_g);
  }

  inline auto &SimulationUnit::get_contribution()
  {
    return this->liquid_scalar->biomass_contribution;
  };

  inline auto &SimulationUnit::getCliq()
  {
    return this->liquid_scalar->C;
  }

  inline auto &SimulationUnit::getCgas()
  {
    return this->gas_scalar->C;
  }

 

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__