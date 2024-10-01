#ifndef __CORE_DATA_EXPORTER_HPP__
#define __CORE_DATA_EXPORTER_HPP__

#include <Kokkos_Core.hpp>
#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

namespace CORE_DE
{

  using export_metadata_t = std::variant<uint64_t, int, std::string>;
  using export_metadata_kv = std::unordered_map<std::string, export_metadata_t>;
  using ViewParticleProperties = Kokkos::View<double **, HostSpace>;
  using simple_export_t =
      std::variant<size_t, std::string, std::vector<size_t>, double>;
  using export_initial_kv = std::unordered_map<std::string, simple_export_t>;

  using matrix_variant_t = std::
      variant<std::span<const double>, std::span<const std::size_t>, double>;

  inline std::string get_filename(std::string_view filename,
                                  std::size_t current_rank)
  {
    return "./results/" + std::string(filename) + "_partial_" +
           std::to_string(current_rank) + ".h5";
  }

  class DataExporter
  {
  public:
    DataExporter(const DataExporter &) = delete;
    DataExporter operator=(const DataExporter &) = delete;

    DataExporter() = default;

  protected:
    struct MultiMatrixDescription
    {
      std::string name;
      std::vector<size_t> dims;
      std::vector<size_t> max_dims;
      std::optional<std::vector<unsigned long long>> chunk_dims;
      bool compression;
      bool is_integer;
    };

    void do_link(std::string_view filename,
                 std::string_view link_name,
                 std::string_view groupname);

    explicit DataExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);
    ~DataExporter();

    void write_matrix(std::string_view name,
                      std::span<const double> values,
                      bool compress = false);

    void write_matrix(std::string_view name,
                      std::span<const double> values,
                      size_t n_row,
                      size_t n_col,
                      bool compress = false);

    void prepare_matrix(MultiMatrixDescription description);

    void append_matrix(std::string_view name, matrix_variant_t data);

    void append_array(std::string_view name,
                      std::span<const double> data,
                      uint64_t last_size = 0);

    void write_properties(std::optional<std::string> specific_dataspace,
                          const export_metadata_kv &values);

    void write_simple(const export_initial_kv &values, std::string_view root);

    void write_simple(std::string specific_dataspace,
                      const simple_export_t &values);

    export_metadata_kv metadata;
    uint64_t export_counter = 0;

  private:
    std::unordered_map<std::string, MultiMatrixDescription> descriptors;
    class impl;
    std::unique_ptr<impl> pimpl;
  };

} // namespace CORE_DE
#endif