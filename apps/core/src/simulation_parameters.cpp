#include <core/simulation_parameters.hpp>
#include <optional>
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
    // params.user_params = UserControlParameters::m_default();

    return params;
  }

  SimulationParameters SimulationParameters::init(const UserControlParameters &user_params)
  {
    // todo
    auto params = SimulationParameters::m_default();
    params.number_exported_result = user_params.number_exported_result;
    params.number_particle = user_params.number_particle;
    params.biomass_initial_concentration = user_params.biomass_initial_concentration;
    params.final_time = user_params.final_time;
    params.results_file_name = user_params.results_file_name;
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
            .serde = false,
            .initialiser_path = "",
            .model_name = "model_light",
            .results_file_name = "",
            .cma_case_path = "",
            .serde_file = std::nullopt};
  }
} // namespace Core