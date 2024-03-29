#ifndef __SYNC_NODES_HPP__
#define __SYNC_NODES_HPP__
#include <simulation/simulation.hpp>
#include <common/execinfo.hpp>

void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation);

void sync_prepare_next(const ExecInfo &exec, Simulation::SimulationUnit &simulation);


#endif //__SYNC_NODES_HPP__