#ifndef __CORE_PARTIAL_EXPORTER_HPP__
#define __CORE_PARTIAL_EXPORTER_HPP__

#include <dataexporter/data_exporter.hpp>
#include <span>

namespace CORE_DE
{
  class PartialExporter : private DataExporter
  {
  public:
    PartialExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);

    void init_fields(uint64_t n_iter, uint64_t n_compartments);

    void write_particle_data(std::span<std::string> names,
                             ViewParticleProperties particle_values,
                             ViewParticleProperties spatial_values,
                             const std::string &ds_name);

    void write_number_particle(const std::vector<size_t> &distribution);

    void write_probe(std::span<const double > data);


  private:
  uint64_t probe_counter_n_element;

  };
}; // namespace CORE_DE

#endif