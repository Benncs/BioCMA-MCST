#ifndef __SIMULATION_PARAMETERS_HPP__
#define __SIMULATION_PARAMETERS_HPP__

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>
namespace Core
{

  /**
   * @struct UserControlParameters
   * @brief A structure to hold user-defined control parameters for simulation settings.
   *
   * This struct is used to define parameters for configuring a simulation. Each
   * parameter allows customization of the simulation's properties, runtime, and output
   * settings.
   */
  struct UserControlParameters
  {

    double biomass_initial_concentration;  ///< Initial concentration of biomass.
    double final_time;                     ///< Final time for the simulation (in seconds).
    double delta_time;                     ///< Time step for the simulation (in seconds).
    uint64_t number_particle;              ///< Number of particles in the simulation.
    int32_t n_thread;                      ///< Number of threads to use for simulation.
    uint32_t number_exported_result;       ///< Number of results to be exported.
    bool recursive;                        ///< Flag to enable recursive processing.
    bool force_override;                   ///< Flag to allow overwriting of existing results.
    bool serde;                            ///< Flag to enable serialization/deserialization.
    std::string initialiser_path;          ///< Path to the initialiser configuration file.
    std::string model_name;                ///< Name of the simulation model.
    std::string results_file_name;         ///< Name of the file where results are saved.
    std::string cma_case_path;             ///< Path to the CMA case file.
    std::optional<std::string> serde_file; ///< Optional file path for serialized data.

    /**
     * @brief Provides default settings for the UserControlParameters structure.
     *
     * This static method returns an instance of UserControlParameters populated
     * with default values for each field.
     */
    static UserControlParameters m_default();
  };

  std::ostream& operator<<(std::ostream& stream, const UserControlParameters& params);

  struct SimulationParameters
  {
    double d_t;
    size_t n_species;
    std::vector<std::string> flow_files;
    uint64_t number_particle;
    bool verbose;
    size_t n_different_maps;
    size_t n_compartments;
    bool is_two_phase_flow = true;
    size_t n_per_flowmap;
    std::string results_file_name;
    double t_per_flow_map;
    uint32_t number_exported_result;
    double biomass_initial_concentration;
    double final_time;
    static SimulationParameters m_default();
    static SimulationParameters init(const UserControlParameters& user_params);
  };

} // namespace Core
#endif //__SIMULATION_PARAMETERS_HPP__
