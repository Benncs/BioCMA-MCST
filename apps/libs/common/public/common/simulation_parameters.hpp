#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
struct SimulationParameters
{
  size_t n_particles;
  size_t n_species;
  double final_time;
  std::vector<std::string> flow_files;
  double d_t = 0;
  int n_threads = -1;
  bool verbose;
  size_t n_different_maps;
  size_t n_compartments;
  bool is_two_phase_flow = true;
  static SimulationParameters m_default();
};

inline std::ostream &operator<<(std::ostream &stream,
                                const SimulationParameters &obj)
{
  stream << obj.n_particles << "\t" << obj.final_time << "\t" << obj.d_t << "\t"
         << obj.n_different_maps << "\t" << obj.n_threads << " ";
  return stream;
}

#endif //__SIMULATION_PARAMETERS_HPP__