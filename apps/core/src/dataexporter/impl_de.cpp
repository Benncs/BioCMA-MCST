

#include <dataexporter/de.hpp>
#include <stdexcept>
#include <string_view>

#define CHECK_PIMPL                                                            \
  if (!pimpl)                                                                  \
  {                                                                            \
    throw std::runtime_error(__FILE__ ": Unexpected ERROR");                   \
  }

#if defined(USE_HIGHFIVE)
#  include <highfive/H5DataSpace.hpp>
#  include <highfive/H5File.hpp>
#  include <highfive/H5PropertyList.hpp>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>

class DataExporter::impl
{
public:
  explicit impl(std::string_view _filename)
      : filename(_filename),
        file(HighFive::File(filename, HighFive::File::Truncate))
  {
  }

  ~impl()
  {
  }

  template <typename T> void write(std::string_view name, T val)
  {
    file.createAttribute<T>(std::string(name), val);
  }

  std::string filename;
  HighFive::File file;
};

DataExporter::DataExporter(std::string_view _filename)
    : pimpl(new impl(_filename))
{
}

void DataExporter::write_properties(
    std::optional<std::string> specific_dataspace,
    const export_metadata_kv &values)
{
  CHECK_PIMPL

  for (const auto &kv : values)
  {
    const std::string &attributeName = kv.first;
    const export_metadata_t &value = kv.second;
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
                                                 // the variant alternative
          pimpl->write<T>(attributeName, val);
        },
        value);
  }
}

void DataExporter::prepare_matrix(
    const std::string &name,
    const std::vector<size_t> &dims,
    const std::vector<size_t> &max_dims,
    const std::vector<unsigned long long> &chunk_dims,
    bool compression)
{
  CHECK_PIMPL

  HighFive::DataSpace dataspace(dims, max_dims);
  HighFive::DataSetCreateProps props;
  props.add(HighFive::Chunking(chunk_dims));
  props.add(HighFive::Shuffle());
  props.add(HighFive::Deflate(compression));
  pimpl->file.createDataSet(name, dataspace, props);
}

#endif