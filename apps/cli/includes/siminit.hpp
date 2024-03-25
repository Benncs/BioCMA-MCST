#ifndef __PROCESS_INIT_HPP__
#define __PROCESS_INIT_HPP__

#include "simulation/models/types.hpp"
#include <common/common.hpp>
#include <flow_iterator.hpp>
#include <simulation/simulation.hpp>
#include <memory>
Simulation::SimulationUnit 
sim_init(ExecInfo &info, SimulationParameters &params, std::shared_ptr<FlowIterator>& it,KModel&& model);

#endif //__PROCESS_INIT_HPP__