#ifndef __HOST_SPECIFIC_HPP__
#define __HOST_SPECIFIC_HPP__

#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <simulation/simulation.hpp>

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  FlowIterator *_flow_handle);

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::shared_ptr<FlowIterator> _flow_handle);

#endif //__HOST_SPECIFIC_HPP__