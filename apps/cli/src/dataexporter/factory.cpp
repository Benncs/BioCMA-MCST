#include <dataexporter/factory.hpp>
#include <dataexporter/data_exporter.hpp>
#include <dataexporter/impl_dataexporter.hpp>


std::unique_ptr<DataExporter>
data_exporter_factory(const ExecInfo &info,
                      const SimulationParameters &params,
                      std::string_view _filename,
                      std::tuple<size_t, size_t> dim,
                      size_t niter,
                      std::span<size_t> distribution)
{
#ifdef USE_HIGHFIVE
  return std::make_unique<DataExportHighFive>(
      info, params, _filename, dim, niter, distribution);
#else
  return std::make_unique<DataExporter>(
      info, params, _filename, dim, niter, distribution);
#endif
}