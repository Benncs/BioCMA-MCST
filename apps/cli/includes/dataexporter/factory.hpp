#ifndef __DATAEXPORT_FACTORY__
#define __DATAEXPORT_FACTORY__

#include <common/simulation_parameters.hpp>
#include <common/execinfo.hpp>
#include <memory>
#include <span>
#include <dataexporter/data_exporter.hpp>

class DataExporter;

std::unique_ptr<DataExporter>
data_exporter_factory(const ExecInfo &info,
                      const SimulationParameters &params,
                      std::string_view _filename,
                      std::tuple<size_t, size_t> dim,
                      size_t niter,
                      std::span<size_t> distribution,double weight=0.);

#endif
