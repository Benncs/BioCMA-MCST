#include <cstddef>
#include <cstdint>
#include <dataexporter/main_exporter.hpp>
#include <optional>
#include <utility>
namespace DE__WIP
{

  MainExporter::MainExporter(const ExecInfo &info,
                             std::string_view _filename,
                             std::optional<export_metadata_t> user_description)
      : DataExporter(info, _filename, std::move(user_description))
  {

    write_properties(std::nullopt, metadata);
    write_simple("misc/n_node_thread", info.thread_per_process);
    write_simple("misc/n_rank", info.n_rank);
  };

  void MainExporter::write_initial(const double weight,
                                   const SimulationParameters &params,
                                   const std::vector<size_t> &distribution)
  {
    export_initial_kv initial_values;
    initial_values["number_particles"] = params.user_params.number_particle;
    initial_values["initial_weight"] = weight;
    initial_values["number_compartment"] = params.n_compartments;
    initial_values["final_time"] = params.user_params.final_time;
    initial_values["particle_distribution"] =
        std::vector<size_t>(distribution.begin(), distribution.end());
    initial_values["delta_time"] = params.d_t;
    initial_values["n_map"] = params.n_different_maps;
    initial_values["t_per_flow_map"] = params.t_per_flow_map;

    write_simple(initial_values);
  }

  void MainExporter::init_fields(uint64_t n_iter,
                                 uint64_t n_compartments,
                                 uint64_t n_species,
                                 bool is_two_phase_flow)
  {
    const uint64_t n_expected_export = n_iter + 2; // Add first + last
    constexpr uint32_t n_growing_arrays_phase = 2; // Concentration + volumes;

    const std::string base_group_name = "records/";

    // First is concentration second is volume
    constexpr std::size_t index_concentration = 0;
    constexpr std::size_t index_volumes = 1;
    auto mK_descriptor = [&](const std::string& phase_name)
    {
      std::vector<unsigned long long> chunk = {1, n_compartments, n_species};
      std::array<MultiMatrixDescription, n_growing_arrays_phase> growing;
      growing[index_concentration] = {
          base_group_name + "concentration_" + phase_name,
          {1, n_compartments, n_species},
          {n_expected_export, n_compartments, n_species},
          chunk,
          true,
          false};

      chunk = {1, n_compartments};

      growing[index_volumes] = {base_group_name + "volume_" + phase_name,
                                {1, n_compartments},
                                {n_expected_export, n_compartments},
                                chunk,
                                true,
                                false};

      return growing;
    };

    for (auto &&descriptor : mK_descriptor("liquid"))
    {
      this->prepare_matrix(descriptor);
    }

    if (is_two_phase_flow)
    {
      for (auto &&descriptor : mK_descriptor("gas"))
      {
        this->prepare_matrix(descriptor);
      }
    }

    this->prepare_matrix(
        {base_group_name + "time", {1}, {n_iter}, std::nullopt, true, false});
  }

  void MainExporter::connect(std::string_view filename,
                             std::string_view groupname)
  {
    do_link(filename, groupname); // TODO
  }

} // namespace DE__WIP