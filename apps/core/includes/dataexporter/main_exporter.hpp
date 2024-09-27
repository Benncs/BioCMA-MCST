#ifndef __CORE_MAIN_EXPORTER_HPP__
#define __CORE_MAIN_EXPORTER_HPP__

#include <common/simulation_parameters.hpp>
#include <dataexporter/ddd_de.hpp>

namespace DE__WIP
{

  class MainExporter : private DataExporter
  {
  public:
    MainExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);

    void write_initial(double weight,
                       const SimulationParameters &params,
                       const std::vector<size_t> &distribution);

    void init_fields(uint64_t n_iter,
                     uint64_t n_compartments,
                     uint64_t n_species,
                     bool is_two_phase_flow);

    void update_fields(double t,
                       std::span<double> concentration_liquid,
                       std::span<const double> liquid_volume,
                       std::optional<std::span<const double>> concentration_gas,
                       std::optional<std::span<const double>> volume_gas);
    void connect(std::string_view filename, std::string_view groupname);
  };

} // namespace DE__WIP

#endif