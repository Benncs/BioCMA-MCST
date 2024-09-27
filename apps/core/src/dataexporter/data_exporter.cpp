#include <Eigen/Core>
#include <chrono>
#include <common/simulation_parameters.hpp>
#include <dataexporter/data_exporter.hpp>
#include <iomanip>
#include <mc/events.hpp>
#include <nl_types.h>
#include <numeric>
#include <simulation/simulation.hpp>
#include <string_view>
#include <tuple>

#ifdef __linux__
#  include <unistd.h>
#endif

static std::string date_time()
{
  std::stringstream ss;
  auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

static std::string get_user_name()
{
  std::string_view res = "someone";
#ifdef __linux__
  char *buff = getlogin();
  if (buff != nullptr)
  {
    res = buff;
  }
#endif
  return std::string(res);
}

DataExporter::DataExporter(const ExecInfo &info,
                           const SimulationParameters &params,
                           std::string_view _filename,
                           std::tuple<size_t, size_t> dim,
                           size_t niter,
                           std::span<size_t> distribution,
                           double weight)
    : filename(_filename), n_row(get<0>(dim)), n_col(get<1>(dim)),
      n_iter(niter + 2)

{
  
  metadata["file_version"] = 3;
  metadata["creation_date"] = date_time();
  metadata["author"] = get_user_name();
  metadata["description"] = "Interesting results";
  metadata["run_id"] = info.run_id;

  initial_values["number_particles"] = params.user_params.number_particle;
  initial_values["initial_weight"] = weight;
  initial_values["number_compartment"] = params.n_compartments;
  initial_values["final_time"] = params.user_params.final_time;
  initial_values["particle_distribution"] =
      std::vector<size_t>(distribution.begin(), distribution.end());
  initial_values["delta_time"] = params.d_t;
  initial_values["n_map"] = params.n_different_maps;
  initial_values["t_per_flow_map"] = params.t_per_flow_map;
}

void DataExporter::write_final_results(Simulation::SimulationUnit &simulation,
                                       std::span<size_t> distribution)
{

  size_t n_part = std::accumulate(
      distribution.begin(), distribution.end(), static_cast<size_t>(0));

  ExportData data = {
      n_part,
      simulation.getCliqData(),
      simulation.getCgasData(),
      &simulation.mc_unit->events,
      simulation.getDim(),
  };

  write_final_results(data, distribution);
}
