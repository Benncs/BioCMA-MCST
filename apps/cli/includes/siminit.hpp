#ifndef __PROCESS_INIT_HPP__
#define __PROCESS_INIT_HPP__

#include "simulation/update_flows.hpp"
#include <memory>

#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <simulation/simulation.hpp>

std::unique_ptr<Simulation::SimulationUnit>
init_simulation(const ExecInfo &info,
                SimulationParameters &params,
                  std::unique_ptr<Simulation::FlowMapTransitioner>& transitioner);

#endif //__PROCESS_INIT_HPP__