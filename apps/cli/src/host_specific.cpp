#include "common/execinfo.hpp"
#include "messages/message_t.hpp"
#include "simulation/simulation.hpp"
#include "simulation/transport.hpp"
#include <common/common.hpp>
#include <host_specific.hpp>
#include <messages/init.hpp>
#include <sync.hpp>
#include <thread>

void computeLiquidFlow(const ExecInfo &info,
                       Simulation::SimulationUnit &unit,
                       FlowInfo &liq_flow)
{
  const auto mat_f_liq =
      FlowmapToMat(liq_flow.flows.data(), liq_flow.flows.getN());
  const auto mat_transition_liq = Simulation::get_transition_matrix(mat_f_liq);
  host_dispatch(info, MPI_W::SIGNALS::RUN, liq_flow.flows.data());

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

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::shared_ptr<FlowIterator> _flow_handle)
{
  ReactorState *f = _flow_handle->get();
  
  double d_t = params.d_t;
  while (f != nullptr)
  {
    simulation.state = f;
    auto &liq_flow = f->liquid_flow;
    auto &gas_flow = f->gas_flow;
    // TODO THREAD POOL
    std::thread _h_liquid(computeLiquidFlow,
                          std::cref(exec),
                          std::ref(simulation),
                          std::ref(liq_flow));
    std::thread _h_gas(
        computeGasFlow, std::ref(simulation), std::ref(gas_flow));

    _h_liquid.join();
    _h_gas.join();

    simulation.cycle_process(d_t);

    sync_step(exec, simulation);
    simulation.step(d_t);
    sync_prepare_next(exec, simulation);
    
    _flow_handle->next(); // this could be done async ?
    f = _flow_handle->get();

    
  }
}