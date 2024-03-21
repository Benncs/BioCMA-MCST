#include "simulation/transport.hpp"
#include <common/common.hpp>
#include <host_specific.hpp>
#include <messages/init.hpp>

#include <iostream>

void computeLiquidFlow(const ExecInfo &info,
                       Simulation::SimulationUnit &unit,
                       FlowInfo &liq_flow)
{
  const auto mat_f_liq =
      FlowmapToMat(liq_flow.flows.data(), liq_flow.flows.getN());
  const auto mat_transition_liq = Simulation::get_transition_matrix(mat_f_liq);
  host_dispatch(info, MPI_SIGNALS::RUN, liq_flow.flows.data());



  unit.setLiquidFlow(
      Simulation::MatFlow(std::move(mat_f_liq), std::move(mat_transition_liq)));
}

void computeGasFlow(Simulation::SimulationUnit &unit, FlowInfo &gas_flow)
{
  const auto mat_f_gas =
      FlowmapToMat(gas_flow.flows.data(), gas_flow.flows.getN());
  const auto mat_transition_gas = Simulation::get_transition_matrix(mat_f_gas);
  unit.setGasFlow(
      Simulation::MatFlow(std::move(mat_f_gas), std::move(mat_transition_gas)));
}