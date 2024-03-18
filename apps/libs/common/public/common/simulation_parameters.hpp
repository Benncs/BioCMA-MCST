#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <string>
#include <vector>

struct SimulationParameters
{
  size_t n_particules;
  size_t n_species;
  double final_time;
  std::vector<std::string> flow_files;
  double d_t = 0;
};

#endif //__SIMULATION_PARAMETERS_HPP__