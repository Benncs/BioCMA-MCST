#ifndef __DATA_EXPORTER_HPP__
#define __DATA_EXPORTER_HPP__

#include "birem_common/macro_constructor_assignment.hpp"
#include "mc/events.hpp"
#include <Eigen/Core>
#include <common/common.hpp>
#include <memory>
#include <simulation/simulation.hpp>
#include <unordered_map>

enum class Dataset
{

};

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
  DataExporter(ExecInfo &info,
               SimulationParameters &params,
               std::string &_filename,
               std::tuple<size_t, size_t> dim,
               size_t niter,
               std::span<size_t> distribution);

  virtual ~DataExporter() = default;
  void write_final_results(const Simulation::SimulationUnit &simulation,
                           std::span<size_t> distribution);

  virtual void append(double t,
                      std::span<double> data,
                      const std::vector<size_t> &distribution,std::span<const double> volume_liquid,std::span<const double> volume_gas){};

  static std::unique_ptr<DataExporter> factory(ExecInfo &info,
                                               SimulationParameters &params,
                                               std::string &_filename,
                                               std::tuple<size_t, size_t> dim,
                                               size_t niter,
                                               std::span<size_t> distribution);

  DELETE_CONSTRUCTORS(DataExporter);
  DELETE_ASSIGMENT(DataExporter)

  size_t n_iter;

protected:
  using export_metadata_t = std::variant<int, std::string>;
  using export_metadata_kv = std::unordered_map<std::string, export_metadata_t>;

  using export_initial_t = std::variant<size_t, std::string,std::vector<size_t>,double>;
  using export_initial_kv = std::unordered_map<std::string, export_initial_t>;

  export_metadata_kv metadata;

  export_initial_kv initial_values;

  std::string filename;
  size_t n_row;
  size_t n_col;
  size_t counter = 0;


  virtual void write_final_results(ExportData &data,
                                   std::span<size_t> distribution)
  {
    std::cerr << "NO implementation specified";
  }
};

#endif //__DATA_EXPORTER_HPP__