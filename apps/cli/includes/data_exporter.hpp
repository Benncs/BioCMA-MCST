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
   DataExporter(ExecInfo& info );
  ~DataExporter() = default;


  DataExporter(DataExporter &&) = delete;
  DataExporter &operator=(DataExporter &&) = delete;
  DataExporter &operator=(const DataExporter &) = delete;
  DataExporter(const DataExporter &) = delete;
};

#endif //__DATA_EXPORTER_HPP__