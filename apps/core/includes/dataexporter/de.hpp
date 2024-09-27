#ifndef __CORE_DATA_EXPORTER_HPP__
#define __CORE_DATA_EXPORTER_HPP__

#include "common/execinfo.hpp"
#include "common/simulation_parameters.hpp"
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>

using export_metadata_t = std::variant<uint64_t, int, std::string>;
using export_metadata_kv = std::unordered_map<std::string, export_metadata_t>;
using ViewParticleProperties = Kokkos::View<double **, HostSpace>;

class DataExporter
{
public:
  explicit DataExporter(std::string_view _filename);
  void write_matrix();
  
  void prepare_matrix(const std::string &name,
                                const std::vector<size_t> &dims,
                                const std::vector<size_t> &max_dims,
                                const std::vector<unsigned long long> &chunk_dims,
                                bool compression);
  void append_matrix();

  void write_properties(std::optional<std::string> specific_dataspace,
                        const export_metadata_kv &values);


  
private:
  class impl; std::unique_ptr<impl> pimpl;
  
};

class MainExporter : public DataExporter
{
public:
  MainExporter(const ExecInfo &info, std::string_view _filename);
  void write_initial(const SimulationParameters &params);
  void update_fields(double t,
                     std::span<double> concentration_liquid,
                     std::span<const double> liquid_volume,
                     std::optional<std::span<const double>> concentration_gas,
                     std::optional<std::span<const double>> volume_gas);
};

class PartialExporter : public DataExporter
{
public:
  PartialExporter(const ExecInfo &info, std::string_view _filename);

  void write_particle_data(std::span<std::string> names,
                           ViewParticleProperties ind_values,
                           ViewParticleProperties spatial_values,
                           const std::string &ds_name);

  void write_number_particle(const std::vector<size_t> &distribution);

private:
};

#endif