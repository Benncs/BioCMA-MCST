#ifndef __DATA_EXPORTER_HPP__
#define __DATA_EXPORTER_HPP__

#include <common/common.hpp>
#include <unordered_map>


enum class Dataset
{
    
};


struct ExportMetaData
{
  
};





class DataExporter
{
public:
   DataExporter(ExecInfo& info ,SimulationParameters& params,std::string_view _filename);
   DataExporter()=delete;
  ~DataExporter() = default;
  void write_final_results();

  DataExporter(DataExporter &&) = delete;
  DataExporter &operator=(DataExporter &&) = delete;
  DataExporter &operator=(const DataExporter &) = delete;
  DataExporter(const DataExporter &) = delete;
  private:
  std::string filename;
};

#endif //__DATA_EXPORTER_HPP__