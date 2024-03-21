#ifndef __HOST_SPECIFIC_HPP__
#define __HOST_SPECIFIC_HPP__

#include "common/execinfo.hpp"
#include "common/simulation_parameters.hpp"
#include "flow_iterator.hpp"
#include "flowmap.hpp"
#include "simulation/simulation.hpp"
#include <simulation/transport.hpp>

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  FlowIterator *_flow_handle);

void computeLiquidFlow(const ExecInfo &info,Simulation::SimulationUnit &unit, FlowInfo &liq_flow);
void computeGasFlow(Simulation::SimulationUnit &unit, FlowInfo &gas_flow);

#endif //__HOST_SPECIFIC_HPP__