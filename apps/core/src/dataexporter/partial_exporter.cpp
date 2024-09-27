#include "dataexporter/ddd_de.hpp"
#include <dataexporter/partial_exporter.hpp>
#include <utility>

namespace DE__WIP
{
  PartialExporter::PartialExporter(
      const ExecInfo &info,
      std::string_view _filename,
      std::optional<export_metadata_t> user_description)
      : DataExporter(info, _filename, std::move(user_description))
  {
    write_properties(std::nullopt, metadata);
  }

  void PartialExporter::init_fields(uint64_t n_iter, uint64_t n_compartments)
  {
    const uint64_t n_expected_export = n_iter + 2; // Add first + last
    MultiMatrixDescription particle_repartition = {
        "records/number_particle",
        {1, n_compartments},
        {n_expected_export, n_compartments},
        std::vector<unsigned long long>({1, n_compartments}),
        true,
        false};

    this->prepare_matrix(particle_repartition);
  }

}; // namespace DE__WIP