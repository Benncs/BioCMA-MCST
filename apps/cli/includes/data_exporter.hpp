#ifndef __DATA_EXPORTER_HPP__
#define __DATA_EXPORTER_HPP__

#include "mc/events.hpp"
#include <common/common.hpp>
#include <unordered_map>
#include <simulation/simulation.hpp>
#include <Eigen/Core>

enum class Dataset
{
    
};


struct ExportMetaData
{
  
};

struct ExportData
{
    size_t number_particles;
    std::span<double> concentration_liquid;
    std::span<double> concentration_gas;
    MC::EventContainer* event;
    std::tuple<size_t,size_t> dim;
};




class DataExporter
{
public:
   DataExporter(ExecInfo& info ,SimulationParameters& params,std::string& _filename,std::tuple<size_t,size_t> dim,size_t niter);
   DataExporter()=delete;
  ~DataExporter() = default;
  void write_final_results(const Simulation::SimulationUnit& simulation,std::span<size_t> distribution);
  void prepare();
  void append(std::span<double> data);

  DataExporter(DataExporter &&) = delete;
  DataExporter &operator=(DataExporter &&) = delete;
  DataExporter &operator=(const DataExporter &) = delete;
  DataExporter(const DataExporter &) = delete;
  private:
  void write_final_results(ExportData& data,std::span<size_t> distribution);
  size_t n_row;
  size_t n_col;
  size_t n_iter;
  std::string filename;
  size_t counter = 0;
};

#endif //__DATA_EXPORTER_HPP__