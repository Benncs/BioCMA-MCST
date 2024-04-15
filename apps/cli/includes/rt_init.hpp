#ifndef __RUNTIME_INIT_HPP__
#define __RUNTIME_INIT_HPP__


#include <common/common.hpp>


void init_environement();



void set_openmp_threads(int rank,
                          int size,
                          ExecInfo &info,
                          SimulationParameters &params);

ExecInfo runtime_init(int argc, char **argv, SimulationParameters &params);

#endif //__RUNTIME_INIT_HPP__

