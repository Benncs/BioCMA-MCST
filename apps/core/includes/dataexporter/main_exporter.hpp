#ifndef __CORE_MAIN_EXPORTER_HPP__
#define __CORE_MAIN_EXPORTER_HPP__

#include "simulation/simulation.hpp"
#include <common/simulation_parameters.hpp>
#include <dataexporter/data_exporter.hpp>

namespace CORE_DE
{

  class MainExporter : private DataExporter
  {
  public:
    MainExporter() = default;
    MainExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);

    void write_initial(double weight,
                       const SimulationParameters &params,
                       const std::vector<size_t> &distribution);

    void write_final(Simulation::SimulationUnit &simulation,
                     std::span<const std::size_t> distribution);

    void init_fields(uint64_t n_iter,
                     uint64_t n_compartments,
                     uint64_t n_species,
                     bool is_two_phase_flow);

    void update_fields(double t,
                       std::span<double> concentration_liquid,
                       std::span<const double> liquid_volume,
                       std::optional<std::span<const double>> concentration_gas,
                       std::optional<std::span<const double>> volume_gas);

    void connect(std::string_view filename,
                 std::string_view link_name,
                 std::string_view groupname);

  private:
    static const std::string base_group_name;
  };

} // namespace CORE_DE

#endif