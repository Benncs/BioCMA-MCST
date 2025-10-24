#ifndef __CORE_PARTIAL_EXPORTER_HPP__
#define __CORE_PARTIAL_EXPORTER_HPP__

#include <common/execinfo.hpp>
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
   * @brief Exporter for each MPI worker (including host) that exports data
   * relative to particle state.
   *
   * This class handles exporting particle-related data, including
   * initialization of fields, writing particle data, managing particle number
   * distributions, and writing probe data. It derives privately from
   * `DataExporter` to encapsulate its backend.
   */
  class PartialExporter : public DataExporter
  {
  public:
    /**
     * @brief Constructs a PartialExporter instance.
     *
     * @param info Execution information necessary for export operations.
     * @param _filename The filename for data export.
     * @param user_description Optional metadata describing the export
     * operation.
     */
    PartialExporter(
        const ExecInfo& info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);

    /**
     * @brief Initializes the fields required for exporting data.
     *
     * @param n_iter Number of iterations for the simulation or process.
     * @param n_compartments Number of compartments to manage.
     */
    void init_fields(uint64_t n_iter, uint64_t n_compartments);

    /**
     * @brief Writes particle data to the output.
     *
     * @param names A span of strings representing particle property names.
     * @param particle_values Particle properties such as mass, age ...
     * @param spatial_values Particle properties such as mass, age ... in space
     * @param ds_name Dataset name associated with this data.
     * @param compress_data Compress particle or not (can be time consuming)
     */
    void write_particle_data(std::span<std::string> names,
                             ViewParticleProperties particle_values,
                             ViewParticleProperties spatial_values,
                             ViewParticleProperties ages_values,
                             const std::string& ds_name,
                             bool compress_data);

    /**
     * @brief Writes the number of particles in each compartment.
     *
     * @param distribution A vector containing particle counts per compartment.
     */
    void write_number_particle(const std::vector<size_t>& distribution);

    /**
     * @brief Writes probe data to the export.
     *
     * @param data A span of constant doubles containing the probe measurements.
     */
    void write_probe(std::span<const double> data);

  private:
    uint64_t probe_counter_n_element; /**< Counter for the number of probe
                                         elements. */
  };

}; // namespace Core

#endif
