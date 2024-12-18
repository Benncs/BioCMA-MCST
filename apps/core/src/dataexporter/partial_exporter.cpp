#include "biocma_cst_config.hpp"
#include "dataexporter/data_exporter.hpp"
#include <dataexporter/partial_exporter.hpp>
#include <utility>

namespace CORE_DE
{
  PartialExporter::PartialExporter(
      const ExecInfo &info,
      std::string_view _filename,
      std::optional<export_metadata_t> user_description)
      : DataExporter(info, _filename, std::move(user_description)),
        probe_counter_n_element(0)
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

    if constexpr (AutoGenerated::FlagCompileTIme::use_probe)
    {

      MultiMatrixDescription probes = {
          "probes",
          {0},
          {(n_expected_export + 1) * AutoGenerated::probe_buffer_size},
          std::vector<unsigned long long>({AutoGenerated::probe_buffer_size}),
          true,
          false};

      this->prepare_matrix(probes);
    }
  }

  void PartialExporter::write_number_particle(
      const std::vector<size_t> &distribution)
  {
    append_matrix("records/number_particle", distribution);

    export_counter++;
  }
  void PartialExporter::write_probe(std::span<const double> data)
  {
    append_array("probes", data, probe_counter_n_element);
    probe_counter_n_element += data.size();
  }

  void
  PartialExporter::write_particle_data(std::span<std::string> names,
                                       ViewParticleProperties particle_values,
                                       ViewParticleProperties spatial_values,
                                       const std::string &ds_name)
  {
    const size_t n_particles = particle_values.extent(1);
    const auto n_compartments = spatial_values.extent(1);
    
    for (size_t i_name = 0; i_name < names.size(); ++i_name)
    {
      auto *ptr_particles =
          Kokkos::subview(particle_values, i_name, Kokkos::ALL).data();

      auto *ptr_spatial =
          Kokkos::subview(spatial_values, i_name, Kokkos::ALL).data();

      this->write_matrix(ds_name + names[i_name], {ptr_particles, n_particles},false);
      this->write_matrix(ds_name + "spatial/" + names[i_name],
                         {ptr_spatial, n_compartments},false);
    }
  }

}; // namespace CORE_DE