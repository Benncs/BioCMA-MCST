#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <cstdint>
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
  std::string model_name;
  double biomass_initial_concentration;
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
  double t_per_flow_map;
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
  stream << obj.user_params.numper_particle << "\t" << obj.user_params.final_time << "\t" << obj.d_t <<
  "\t"
         << obj.n_different_maps << "\t" << obj.user_params.n_thread << " ";
  return stream;
}

#endif //__SIMULATION_PARAMETERS_HPP__