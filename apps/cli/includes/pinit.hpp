#ifndef __PROCESS_INIT_HPP__
#define __PROCESS_INIT_HPP__

#include <common/common.hpp>
#include <flow_iterator.hpp>
#include <simulation/simulation.hpp>

SimulationUnit
pinit(ExecInfo &info, SimulationParameters &params, FlowIterator **it);

#endif //__PROCESS_INIT_HPP__