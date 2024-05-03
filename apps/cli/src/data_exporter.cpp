#include "common/simulation_parameters.hpp"
#include <data_exporter.hpp>

#include <chrono>

#include <iomanip>

std::string date_time()
{
  std::stringstream ss;
  auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}
#ifdef USE_HIGHFIVE
#  include <Eigen/Dense>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>
static void write_attributes(HighFive::File &file)
{

  file.createAttribute<int>("file_version", 1);
  file.createAttribute<std::string>("creation_date", date_time());
  file.createAttribute<std::string>("author", std::string("someone"));
  file.createAttribute<std::string>("description",
                                    std::string("Interesting results"));
}

static void write_initial(HighFive::File &file,
                          ExecInfo &info,
                          SimulationParameters &params)
{
  size_t n_0_particle = info.n_rank * params.n_particles;

  file.createDataSet("initial_parameters/number_particles_0", n_0_particle);
  file.createDataSet("initial_parameters/number_compartment",
                     params.n_compartments);
  file.createDataSet("initial_parameters/final_time", params.final_time);
  file.createDataSet("initial_parameters/delta_time", params.d_t);
  file.createDataSet("initial_parameters/n_map", params.n_different_maps);
  file.createDataSet("misc/n_node_thread", info.thread_per_process);
  file.createDataSet("misc/n_rank", info.n_rank);
}

DataExporter::DataExporter(ExecInfo &info,
                           SimulationParameters &params,
                           std::string_view _filename)
    : filename(_filename.data())

{
  HighFive::File file(filename, HighFive::File::Truncate);
  write_initial(file, info, params);
  write_attributes(file);
}

#else
DataExporter::DataExporter(ExecInfo &info,
                           SimulationParameters &params,
                           std::string_view _filename)
    : filename(_filename.data())

{
  // NOP
}
#endif