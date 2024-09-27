#ifndef __CORE_PARTIAL_EXPORTER_HPP__
#define __CORE_PARTIAL_EXPORTER_HPP__

#include <dataexporter/ddd_de.hpp>
#include <span>

namespace DE__WIP
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
                             ViewParticleProperties ind_values,
                             ViewParticleProperties spatial_values,
                             const std::string &ds_name);

    void write_number_particle(const std::vector<size_t> &distribution);

  private:
  };
}; // namespace DE__WIP

#endif