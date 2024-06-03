#ifndef __HOST_SPECIFIC_HPP__
#define __HOST_SPECIFIC_HPP__

#include "data_exporter.hpp"
#include "simulation/update_flows.hpp"
#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <simulation/simulation.hpp>

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  std::unique_ptr<Simulation::FlowMapTransitioner>&& transitioner);

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,DataExporter* exporter);

#endif //__HOST_SPECIFIC_HPP__