#ifndef __CORE_DATA_EXPORTER_HPP__
#define __CORE_DATA_EXPORTER_HPP__

#include <Kokkos_Core.hpp>
#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

namespace DE__WIP
{

  using export_metadata_t = std::variant<uint64_t, int, std::string>;
  using export_metadata_kv = std::unordered_map<std::string, export_metadata_t>;
  using ViewParticleProperties = Kokkos::View<double **, HostSpace>;
  using simple_export_t =
      std::variant<size_t, std::string, std::vector<size_t>, double>;
  using export_initial_kv = std::unordered_map<std::string, simple_export_t>;
  class DataExporter
  {
  public:
    DataExporter(const DataExporter &) = delete;
    DataExporter(DataExporter &&) = delete;
    DataExporter operator=(DataExporter &&) = delete;
    DataExporter operator=(const DataExporter &) = delete;

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

    void do_link(std::string_view filename,std::string_view groupname);

    explicit DataExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);
    ~DataExporter();

    void write_matrix();

    void prepare_matrix(MultiMatrixDescription description);

    void append_matrix();

    void write_properties(std::optional<std::string> specific_dataspace,
                          const export_metadata_kv &values);

    void write_simple(const export_initial_kv &values);

    void write_simple(std::string specific_dataspace,
                      const simple_export_t &values);

    export_metadata_kv metadata;

  private:
    class impl;
    std::unique_ptr<impl> pimpl;
  };

} // namespace DE__WIP
#endif