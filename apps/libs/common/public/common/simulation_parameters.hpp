#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>



struct UserControlParameters
{
  std::string cma_case_path;
  bool recursive;
  double final_time;
  uint64_t numper_particle;
  double delta_time;
  std::string results_file_name;
  int n_thread;
  uint32_t number_exported_result;
  static UserControlParameters m_default();
};



struct SimulationParameters
{

  UserControlParameters user_params;
  double d_t;
  size_t n_species;
  std::vector<std::string> flow_files;
  bool verbose;
  size_t n_different_maps;
  size_t n_compartments;
  bool is_two_phase_flow = true;
  size_t n_per_flowmap;
  std::string results_file_name;
  static SimulationParameters m_default();
};

struct ExportParameters
{
  size_t n_save; 
  std::string filename; 
};

inline std::ostream &operator<<(std::ostream &stream,
                                const SimulationParameters &obj)
{
  // stream << obj.n_particles << "\t" << obj.final_time << "\t" << obj.d_t << "\t"
  //        << obj.n_different_maps << "\t" << obj.n_threads << " ";
  return stream;
}

#endif //__SIMULATION_PARAMETERS_HPP__