// DEFAULT
#ifndef USE_HIGHFIVE
#  include <dataexporter/data_exporter.hpp>

namespace Core
{
  class DataExporter::impl
  {
  };
  void DataExporter::do_link(std::string_view filename,
                             std::string_view link_name,
                             std::string_view groupname)
  {
  }

  DataExporter::DataExporter(const ExecInfo& info,
                             std::string_view _filename,
                             std::optional<export_metadata_t> user_description)
  {
    if (logger)
    {
      logger->error("NOT DATAEXPORTER IMPLEMENTATION");
    }
  }
  DataExporter::~DataExporter() = default;

  void DataExporter::write_matrix(std::string_view name,
                                  std::span<const double> values,
                                  bool compress)
  {
  }

  void DataExporter::prepare_matrix(MultiMatrixDescription description)
  {
  }

  void DataExporter::append_matrix(std::string_view name, matrix_variant_t data)
  {
  }

  void DataExporter::append_array(std::string_view name,
                                  std::span<const double> data,
                                  uint64_t last_size)
  {
  }

  void
  DataExporter::write_properties(std::optional<std::string> specific_dataspace,
                                 const export_metadata_kv& values)
  {
  }

  void DataExporter::write_matrix(std::string_view name,
                                  std::span<const double> values,
                                  size_t n_row,
                                  size_t n_col,
                                  bool compress)
  {
  }

  void DataExporter::write_simple(const export_initial_kv& values,
                                  std::string_view root)
  {
  }

  void DataExporter::write_simple(std::string specific_dataspace,
                                  const simple_export_t& values)
  {
  }
} // namespace Core
#endif
