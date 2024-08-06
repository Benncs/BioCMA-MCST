
#ifdef USE_HIGHFIVE
#  include "dataexporter/data_exporter.hpp"
#  ifndef __IMPL_DATAEXPORTER_HPP__
#    define __IMPL_DATAEXPORTER_HPP__

namespace HighFive
{
  class File;
} // namespace HighFive

class DataExportHighFive : public DataExporter
{
public:
  DataExportHighFive(const ExecInfo &info,
                     const SimulationParameters &params,
                     std::string_view _filename,
                     std::tuple<size_t, size_t> dim,
                     size_t niter,
                     std::span<size_t> distribution);

  SET_NON_COPYABLE(DataExportHighFive)
  SET_NON_MOVABLE(DataExportHighFive)

  ~DataExportHighFive()
  {
    // delete _file;
  }

  void append(double t,
              std::span<double> concentration_liquid,
              const std::vector<size_t> &distribution,
              std::span<const double> liquid_volume,
              std::span<const double> volume_gas) override;

  void write_final_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>> &,
      const std::unordered_map<std::string, std::vector<double>> &) override;

  void write_initial_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>>
          & /*unused*/,
      const std::unordered_map<std::string, std::vector<double>> & /*unused*/)
      override;



protected:
  void write_final_results(ExportData &data,
                           std::span<size_t> distribution) override;

private:
  void prepare();
  static void write_attributes(HighFive::File &file, export_metadata_kv &md);

  static void write_initial(HighFive::File &file,
                            const ExecInfo &info,
                            export_initial_kv &md);

  void write_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>> &,
      const std::unordered_map<std::string, std::vector<double>> &,
      const std::string &ds_name);

  // HighFive::File *_file;
};

#  endif

#endif