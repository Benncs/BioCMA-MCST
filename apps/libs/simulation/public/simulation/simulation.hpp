#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include "common/execinfo.hpp"
#include <common/common.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <random>
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

    std::unique_ptr<MC::MonteCarloUnit> unit;
    std::unique_ptr<MC::ParticlesContainer> container;
    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;

    void post_init(KModel &&_km);

    inline auto &getC()
    {
      return this->liquid_scalar->C;
    }

    void setVolumeLiquid(std::vector<double> &&volumes)
    {
      //   this->liquid_scalar->V = Eigen::Map<Eigen::DiagonalMatrix<double,
      //   Eigen::Dynamic>>(
      // volumes.data(), this->unit->domain.n_compartments(),
      // this->unit->domain.n_compartments());
      this->liquid_scalar->V.diagonal() =
          Eigen::Map<const Eigen::VectorXd>(volumes.data(), volumes.size());
    }

    void step(double d_t);
    void cycle_process(const double d_t);

    inline void setLiquidFlow(MatFlow &&_flows_l)
    {
      flow_liquid = std::move(_flows_l);
    }

    inline void setGasFlow(MatFlow &&_flows_g)
    {
      flow_gas = std::move(_flows_g);
    }

    auto &get_contribution()
    {
      return this->liquid_scalar->biomass_contribution;
    };

    void clear_contribution()
    {
      for (auto &&i : contribs)
      {
        i.setZero();
      }
      this->liquid_scalar->biomass_contribution.setZero();

      for (auto i = 0; i < liquid_scalar->C.cols(); ++i)
      {
        liquid_scalar->C.coeffRef(1, i) = 8e-3;
      }
    }

    MatFlow flow_liquid; // TODO OPTI
    MatFlow flow_gas;    // TODO OPTI
  private:
    void post_init_container();
    void post_init_compartments();
    bool host;
    KModel kmodel;
    size_t np;
    std::unique_ptr<ScalarSimulation> liquid_scalar;
    std::unique_ptr<ScalarSimulation> gas_scalar;

    std::vector<Eigen::MatrixXd> contribs;
  };
} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__