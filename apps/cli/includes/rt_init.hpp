#ifndef __RUNTIME_INIT_HPP__
#define __RUNTIME_INIT_HPP__

#include <common/common.hpp>
#include <data_exporter.hpp>

void init_environment();
std::string env_file_path();
void append_date_time(std::ofstream &fd);

void set_openmp_threads(int rank,
                        int size,
                        ExecInfo &info,
                        SimulationParameters &params);

ExecInfo runtime_init(int argc, char **argv, SimulationParameters &params);

#endif //__RUNTIME_INIT_HPP__
