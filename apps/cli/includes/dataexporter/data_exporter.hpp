#ifndef __DATA_EXPORTER_HPP__
#define __DATA_EXPORTER_HPP__

#include <cmt_common/macro_constructor_assignment.hpp>
#include <common/common.hpp>
#include <mc/events.hpp>
#include <mc/particles/particle_model.hpp>
#include <string_view>
#include <unordered_map>

// Foraward declaration
namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

struct ExportData
{
  size_t number_particles;
  std::span<double> concentration_liquid;
  std::span<double> concentration_gas;
  MC::EventContainer *event;
  std::tuple<size_t, size_t> dim;
};

class DataExporter
{
public:
  DataExporter(const ExecInfo &info,
               const SimulationParameters &params,
               std::string_view _filename,
               std::tuple<size_t, size_t> dim,
               size_t niter,
               std::span<size_t> distribution,
               double weight);

  virtual ~DataExporter() = default;
  void write_final_results(Simulation::SimulationUnit &simulation,
                           std::span<size_t> distribution);

  virtual void write_final_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>>
          &props,
      const std::unordered_map<std::string, std::vector<double>>
          &spatial_props) = 0;

  virtual void write_initial_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>>
          &props,
      const std::unordered_map<std::string, std::vector<double>>
          &spatial_props) = 0;

  virtual void append(double t,
                      std::span<double> concentration_liquid,
                      const std::vector<size_t> &distribution,
                      std::span<const double> liquid_volume,
                      std::span<const double> volume_gas) = 0;

  virtual void append_particle_properties(
      size_t counter,
      const std::unordered_map<std::string, std::vector<model_properties_t>>
          &props,
      const std::unordered_map<std::string, std::vector<double>>
          &spatial_props) = 0;

  DELETE_CONSTRUCTORS(DataExporter)
  DELETE_ASSIGMENT(DataExporter)
  size_t expectecNiteration() const
  {
    return n_iter;
  }

protected:
  using export_metadata_t = std::variant<int, std::string>;
  using export_metadata_kv = std::unordered_map<std::string, export_metadata_t>;

  using export_initial_t =
      std::variant<size_t, std::string, std::vector<size_t>, double>;
  using export_initial_kv = std::unordered_map<std::string, export_initial_t>;

  export_metadata_kv metadata;

  export_initial_kv initial_values;

  std::string filename;
  size_t n_row;
  size_t n_col;
  size_t counter = 0;
  size_t n_iter;

  virtual void write_final_results(ExportData & /*data*/,
                                   std::span<size_t> /*distribution*/)
  {
    std::cerr << "NO implementation specified";
  }
};

#endif //__DATA_EXPORTER_HPP__