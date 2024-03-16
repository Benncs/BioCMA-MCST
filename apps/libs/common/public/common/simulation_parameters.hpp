#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <string>
#include <vector>

struct SimulationParameters
{
  size_t n_particules;
  double final_time;
  std::vector<std::string> flow_files;
};

#endif //__SIMULATION_PARAMETERS_HPP__