#ifndef __PROCESS_INIT_HPP__
#define __PROCESS_INIT_HPP__

#include <common/common.hpp>
#include <cma_read/flow_iterator.hpp>
#include <simulation/simulation.hpp>
#include <memory>

Simulation::SimulationUnit
init_simulation(ExecInfo &info,
                SimulationParameters &params,
                std::shared_ptr<FlowIterator> &_flow_handle,
                KModel &&model);

#endif //__PROCESS_INIT_HPP__