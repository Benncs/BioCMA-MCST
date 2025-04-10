#include <core/simulation_parameters.hpp>
#include <optional>
#include <ostream>
namespace Core
{

  SimulationParameters SimulationParameters::m_default()
  {
    SimulationParameters params;
    params.n_species = 1;
    params.flow_files = {"./cma_data/"};
    params.verbose = false;
    params.n_compartments = 0;
    params.t_per_flow_map = 0;
    params.save_final_serde=false;
    // params.user_params = UserControlParameters::m_default();

    return params;
  }

  SimulationParameters SimulationParameters::init(const UserControlParameters& user_params)
  {
    // todo
    
    auto params = SimulationParameters::m_default();
    params.number_exported_result = user_params.number_exported_result;
    params.number_particle = user_params.number_particle;
    params.biomass_initial_concentration = user_params.biomass_initial_concentration;
    params.final_time = user_params.final_time;
    params.results_file_name = user_params.results_file_name;
    params.save_final_serde =user_params.save_serde;
    return params;
  }

  UserControlParameters UserControlParameters::m_default()
  {
    constexpr double x0 = 1.;
    return {.biomass_initial_concentration = x0,
            .final_time = 0,
            .delta_time = 0.,
            .number_particle = 1,
            .n_thread = -1,
            .number_exported_result = 0,
            .recursive = false,
            .force_override = false,
            .load_serde = false,
            .save_serde=false,
            .initialiser_path = "",
            .model_name = "None",
            .results_file_name = "",
            .cma_case_path = "",
            .serde_file = std::nullopt,
            };
  }
  std::ostream& operator<<(std::ostream& stream, const UserControlParameters& params)
  {
    stream << "UserControlParameters:\n"
           << "  Biomass Initial Concentration: " << params.biomass_initial_concentration << "\n"
           << "  Final Time: " << params.final_time << "\n"
           << "  Delta Time: " << params.delta_time << "\n"
           << "  Number of Particles: " << params.number_particle << "\n"
           << "  Number of Threads: " << params.n_thread << "\n"
           << "  Number of Exported Results: " << params.number_exported_result << "\n"
           << "  Recursive: " << (params.recursive ? "true" : "false") << "\n"
           << "  Force Override: " << (params.force_override ? "true" : "false") << "\n"
           << "  Serde: " << (params.load_serde ? "true" : "false") << "\n"
           << "  Initialiser Path: " << params.initialiser_path << "\n"
           << "  Model Name: " << params.model_name << "\n"
           << "  Results File Name: " << params.results_file_name << "\n"
           << "  CMA Case Path: " << params.cma_case_path << "\n"
           << "  Serde File: "
           << (params.serde_file.has_value() ? params.serde_file.value() : "nullopt") << "\n";
    return stream;
  }

} // namespace Core
