#ifndef __RUNTIME_INIT_HPP__
#define __RUNTIME_INIT_HPP__

#include <common/common.hpp>
#include <data_exporter.hpp>
#include <chrono> 
#include <iomanip>

void init_environment();
std::string env_file_path();
// void append_date_time(std::ofstream &fd);


// std::string append_date_time(std::string_view string);



template<typename Stream>
void append_date_time(Stream &stream)
{
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    stream << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
}

std::string sappend_date_time(std::string_view string);


void set_openmp_threads(int rank,
                        int size,
                        ExecInfo &info,
                        const UserControlParameters &params);

ExecInfo runtime_init(int argc, char **argv, const SimulationParameters &params);

void register_run(const ExecInfo &exec, SimulationParameters &params);

#endif //__RUNTIME_INIT_HPP__
