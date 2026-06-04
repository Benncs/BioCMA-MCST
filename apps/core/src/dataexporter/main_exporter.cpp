#include <biocma_cst_config.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
#include <cstdint>
#include <dataexporter/main_exporter.hpp>
#include <filesystem>
#include <iomanip>
#include <optional>
#include <simulation/descriptors/dimensions.hpp>
#include <simulation/simulation_getter.hpp>
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

  void
  MainExporter::write_initial(const double weight,
                              const Core::SimulationParameters& params)
  {
    export_initial_kv initial_values;
    initial_values["number_particles"] = params.number_particle;
    initial_values["initial_weight"] = weight;
    initial_values["initial_biomass_concentration"]
        = params.biomass_initial_concentration;
    initial_values["number_compartment"] = params.n_compartments;
    initial_values["final_time"] = params.final_time;
    // initial_values["particle_distribution"] =
    //     std::vector<size_t>(distribution.begin(), distribution.end());
    initial_values["delta_time"] = params.d_t;
    initial_values["n_map"] = params.n_different_maps;
    initial_values["t_per_flow_map"] = params.t_per_flow_map;

    write_simple(initial_values, "initial_parameters/");
  }

  void
  MainExporter::init_fields(uint64_t n_iter,
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
      std::vector<unsigned long long> chunk = { 1, n_compartments, n_species };
      std::array<MultiMatrixDescription, n_growing_arrays_phase> growing;

      growing[index_concentration]
          = { .name = base_group_name + "concentration_" + phase_name,
              .dims = { 1, n_compartments, n_species },
              .max_dims = { n_expected_export, n_compartments, n_species },
              .chunk_dims = chunk,
              .compression = true,
              .is_integer = false };

      chunk = { 1, n_compartments };

      growing[index_volumes]
          = { .name = base_group_name + "volume_" + phase_name,
              .dims = { 1, n_compartments },
              .max_dims = { n_expected_export, n_compartments },
              .chunk_dims = chunk,
              .compression = true,
              .is_integer = false };

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

      std::vector<unsigned long long> chunk = { 1, n_compartments, n_species };
      this->prepare_matrix({ base_group_name + "mtr",
                             { 1, n_compartments, n_species },
                             { n_expected_export, n_compartments, n_species },
                             chunk,
                             true,
                             false });
    }

    this->prepare_matrix({ base_group_name + "time",
                           { 1 },
                           { n_expected_export },
                           std::vector<unsigned long long>({ 1 }),
                           true,
                           false });
  }

  void
  MainExporter::update_fields(
      double t,
      std::span<const double> concentration_liquid,
      std::span<const double> liquid_volume,
      std::optional<std::span<const double>> concentration_gas,
      std::optional<std::span<const double>> volume_gas,
      std::optional<std::span<const double>> mtr)
  {

    if (export_counter >= n_expected_export)
    {
      return;
    }
    append_matrix(base_group_name + "concentration_liquid",
                  concentration_liquid);
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

  void
  MainExporter::write_final(const Simulation::Getter& getter,
                            std::size_t number_particles)
  {
    const auto& event = getter.mc_unit()->events;

    export_initial_kv final_values;
    final_values["number_particles"] = number_particles;

    if constexpr (AutoGenerated::FlagCompileTime::enable_event_counter)
    {
      final_values["events/move"] = event.get<MC::EventType::Move>();
      final_values["events/total_division"]
          = event.get<MC::EventType::NewParticle>();
      final_values["events/total_death"] = event.get<MC::EventType::Death>();
      final_values["events/total_exit"] = event.get<MC::EventType::Exit>();
    }

    const auto [n_row, n_col] = getter.getDimensions();

    write_simple(final_values, "final_result/");

    // Compress data is encouraged, as final is only written once let always
    // compress
    const bool f_compress = true;

    write_matrix("final_result/concentration_liquid",
                 getter.getCliqData(),
                 n_row,
                 n_col,
                 f_compress);

    auto opt_gas = getter.getCgasData();
    if (opt_gas.has_value())
    {
      write_matrix(
          "final_result/concentration_gas", *opt_gas, n_row, n_col, f_compress);
    }
  }

  bool
  fill_and_check_result_file_path(const std::shared_ptr<IO::Logger>& logger,
                                  Core::UserControlParameters& params)
  {

    auto sappend_date_time = [](std::string_view string) noexcept
    {
      std::stringstream fd;
      fd << string;
      auto now = std::chrono::system_clock::to_time_t(
          std::chrono::system_clock::now());
      fd << std::put_time(std::localtime(&now), "%Y-%m-%d-%H:%M:%S");
      return fd.str();
    };
    // TODO
    try
    {

      std::filesystem::path path;

      if (params.results_file_name.empty())
      {
        // TODO
        const auto time_stamped_folder_name = sappend_date_time("result_");
        path = std::filesystem::canonical(AutoGenerated::result_dir_path);
        path.append(time_stamped_folder_name);
        // In this specific case, create directory
        if (logger)
        {
          logger->alert(
              "Deprecated",
              "Result path not given, using default result path will be "
              "deprecated soon");
        }

        path.append(time_stamped_folder_name);
      }
      else
      {
        path.append(params.results_file_name);
      }
      path = std::filesystem::weakly_canonical(path);
      params.results_file_name = path;

      return std::filesystem::exists(path.parent_path());
    }
    catch (std::exception& e)
    {
      if (logger)
      {
        logger->error(
            "Result path not given, using default result path will be "
            "deprecated soon");
      }

      return false;
    }
  }

} // namespace Core
