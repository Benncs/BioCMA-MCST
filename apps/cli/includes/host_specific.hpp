#ifndef __HOST_SPECIFIC_HPP__
#define __HOST_SPECIFIC_HPP__

#include <dataexporter/factory.hpp>
#include "simulation/update_flows.hpp"
#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <memory>
#include <simulation/simulation.hpp>

void host_process(const ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  const SimulationParameters &params,
                  std::unique_ptr<Simulation::FlowMapTransitioner>&& transitioner);

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,std::unique_ptr<DataExporter>& exporter);

#endif //__HOST_SPECIFIC_HPP__