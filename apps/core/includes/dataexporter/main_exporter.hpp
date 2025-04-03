#ifndef __CORE_MAIN_EXPORTER_HPP__
#define __CORE_MAIN_EXPORTER_HPP__

#include "mc/events.hpp"
#include "simulation/simulation.hpp"
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdint>
#include <dataexporter/data_exporter.hpp>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>
namespace Core
{
  /**
   * @brief Exporter for host-only data related to general-purpose simulation data
   * (e.g., concentration, parameters, and other global data).
   *
   * This class is responsible for exporting simulation data at both the initial
   * and final stages, as well as for managing field updates during the simulation.
   * It provides functionalities to handle both single-phase and two-phase flow simulations.
   */
  class MainExporter : public DataExporter
  {
  public:
    /**
     * @brief Default constructor for MainExporter.
     */
    MainExporter() = default;

    /**
     * @brief Constructs a MainExporter instance.
     *
     * @param info Execution information required for export operations.
     * @param _filename The filename for data export.
     * @param user_description Optional metadata describing the export.
     */
    MainExporter(const ExecInfo& info,
                 std::string_view _filename,
                 std::optional<export_metadata_t> user_description = std::nullopt);

    /**
     * @brief Writes initial simulation data to the output.
     *
     * @param weight Initial particle weight used in the simulation.
     * @param params Core simulation parameters
     * @param distribution Initial distribution of particles or entities across compartments.
     */
    void write_initial(double weight,
                       const Core::SimulationParameters& params);

    /**
     * @brief Writes final simulation data to the output.
     *
     * @param simulation Reference to the simulation unit containing final state information.
     * @param distribution Final distribution of particles or entities across compartments.
     */
    void write_final(Simulation::SimulationUnit& simulation,std::size_t number_particles);

    /**
     * @brief Initializes fields required for simulation export.
     *
     * @param n_iter Number of iterations for the simulation.
     * @param n_compartments Number of compartments to manage.
     * @param n_species Number of species or particle types being simulated.
     * @param is_two_phase_flow Flag indicating if the simulation involves two-phase flow.
     */
    void init_fields(uint64_t n_iter,
                     uint64_t n_compartments,
                     uint64_t n_species,
                     bool is_two_phase_flow);

    /**
     * @brief Updates fields with the current simulation state.
     *
     * @param t Current simulation time.
     * @param concentration_liquid A span of doubles representing liquid phase concentrations.
     * @param liquid_volume A span of doubles representing liquid phase volumes.
     * @param concentration_gas Optional span of doubles for gas phase concentrations.
     * @param volume_gas Optional span of doubles for gas phase volumes.
     */
    void update_fields(double t,
                       std::span<double> concentration_liquid,
                       std::span<const double> liquid_volume,
                       std::optional<std::span<const double>> concentration_gas,
                       std::optional<std::span<const double>> volume_gas,
                       std::optional<std::span<const double>> mtr,std::optional<std::span<std::size_t>> events);

  private:
    static const std::string
        base_group_name;          /**< Base group name for exported data organization. */
    uint64_t n_expected_export{}; /**< Expected number of exports for tracking purposes. */
  };

} // namespace Core

#endif