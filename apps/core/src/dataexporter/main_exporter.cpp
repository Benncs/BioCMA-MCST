#include "biocma_cst_config.hpp"
#include <common/execinfo.hpp>
#include <cstddef>
#include <cstdint>
#include <dataexporter/main_exporter.hpp>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace Core
{
  const std::string MainExporter::base_group_name = "records/";
  MainExporter::MainExporter(const ExecInfo& info,
                             std::string_view _filename,
                             std::optional<export_metadata_t> user_description)
      : DataExporter(info, _filename, std::move(user_description))
  {

    write_properties(std::nullopt, metadata);
    write_simple("misc/n_node_thread", info.thread_per_process);
    write_simple("misc/n_rank", info.n_rank);
  };

  void MainExporter::write_initial(const double weight, const Core::SimulationParameters& params)
  {
    export_initial_kv initial_values;
    initial_values["number_particles"] = params.number_particle;
    initial_values["initial_weight"] = weight;
    initial_values["initial_biomass_concentration"] = params.biomass_initial_concentration;
    initial_values["number_compartment"] = params.n_compartments;
    initial_values["final_time"] = params.final_time;
    // initial_values["particle_distribution"] =
    //     std::vector<size_t>(distribution.begin(), distribution.end());
    initial_values["delta_time"] = params.d_t;
    initial_values["n_map"] = params.n_different_maps;
    initial_values["t_per_flow_map"] = params.t_per_flow_map;

    write_simple(initial_values, "initial_parameters/");
  }

  void MainExporter::init_fields(uint64_t n_iter,
                                 uint64_t n_compartments,
                                 uint64_t n_species,
                                 bool is_two_phase_flow)
  {
    n_expected_export = n_iter + 2;                // Add first + last
    constexpr uint32_t n_growing_arrays_phase = 2; // Concentration + volumes;
    // First is concentration second is volume
    constexpr std::size_t index_concentration = 0;
    constexpr std::size_t index_volumes = 1;
    auto mK_descriptor = [&](const std::string& phase_name)
    {
      std::vector<unsigned long long> chunk = {1, n_compartments, n_species};
      std::array<MultiMatrixDescription, n_growing_arrays_phase> growing;
      growing[index_concentration] = {base_group_name + "concentration_" + phase_name,
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

    for (auto&& descriptor : mK_descriptor("liquid"))
    {
      this->prepare_matrix(descriptor);
    }

    if (is_two_phase_flow)
    {
      for (auto&& descriptor : mK_descriptor("gas"))
      {
        this->prepare_matrix(descriptor);
      }

      std::vector<unsigned long long> chunk = {1, n_compartments, n_species};
      this->prepare_matrix({base_group_name + "mtr",
                            {1, n_compartments, n_species},
                            {n_expected_export, n_compartments, n_species},
                            chunk,
                            true,
                            false});
    }

    this->prepare_matrix({base_group_name + "time",
                          {1},
                          {n_expected_export},
                          std::vector<unsigned long long>({1}),
                          true,
                          false});
  }

  void MainExporter::update_fields(double t,
                                   std::span<double> concentration_liquid,
                                   std::span<const double> liquid_volume,
                                   std::optional<std::span<const double>> concentration_gas,
                                   std::optional<std::span<const double>> volume_gas,
                                   std::optional<std::span<const double>> mtr)
  {
    if (export_counter >= n_expected_export)
    {
      return;
    }
    append_matrix(base_group_name + "concentration_liquid", concentration_liquid);
    append_matrix(base_group_name + "volume_liquid", liquid_volume);
    append_matrix(base_group_name + "time", t);
    if (concentration_gas.has_value())
    {
      append_matrix(base_group_name + "concentration_gas", *concentration_gas);
    }
    if (volume_gas.has_value())
    {
      append_matrix(base_group_name + "volume_gas", *volume_gas);
    }
    if (mtr.has_value())
    {
      append_matrix(base_group_name + "mtr", *mtr);
    }

    export_counter++;
  }

  void MainExporter::write_final(Simulation::SimulationUnit& simulation,
                                 std::size_t number_particles)
  {
    ;
    const auto& event = simulation.mc_unit->events;

    export_initial_kv final_values;
    final_values["number_particles"] = number_particles;
    // final_values["distribution"] = std::vector<std::size_t>(
    //     distribution.begin(), distribution.end()); // TODO Remove ?
    if constexpr (AutoGenerated::FlagCompileTime::enable_event_counter)
    {
      final_values["events/move"] = event.get<MC::EventType::Move>();
      final_values["events/total_division"] = event.get<MC::EventType::NewParticle>();
      final_values["events/total_death"] = event.get<MC::EventType::Death>();
      final_values["events/total_exit"] = event.get<MC::EventType::Exit>();
    }

    const auto [n_row, n_col] = simulation.getDimensions();

    write_simple(final_values, "final_result/");

    write_matrix("final_result/concentration_liquid", simulation.getCliqData(), n_row, n_col, true);

    auto opt_gas = simulation.getCgasData();
    if (opt_gas.has_value())
    {
      write_matrix("final_result/concentration_gas", *opt_gas, n_row, n_col, true);
    }
  }

} // namespace Core