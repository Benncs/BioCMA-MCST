#ifndef __RUNTIME_INIT_HPP__
#define __RUNTIME_INIT_HPP__


#include "common/simulation_parameters.hpp"
#include <common/execinfo.hpp>

void set_openmp_threads(int rank,
                          int size,
                          ExecInfo &info,
                          SimulationParameters &params);

ExecInfo runtime_init(int argc, char **argv, SimulationParameters &params);

#endif //__RUNTIME_INIT_HPP__

