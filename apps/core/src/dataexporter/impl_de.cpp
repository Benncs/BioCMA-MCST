#include <dataexporter/ddd_de.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/eigen.hpp>
#include <highfive/highfive.hpp>
#include <stdexcept>
#include <string_view>

#ifdef __linux__
#  include <pwd.h>
#  include <unistd.h>
#endif

#define CHECK_PIMPL                                                            \
  if (!pimpl)                                                                  \
  {                                                                            \
    throw std::runtime_error(__FILE__ ": Unexpected ERROR");                   \
  }
namespace DE__WIP
{
  static std::string date_time()
  {
    // Non vedo l’ora che arrivi c++23-format
    std::stringstream ss;
    auto now =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
    return ss.str();
  }

  static std::string get_user_name()
  {
    std::string_view res = "someone";
#ifdef __linux__
    uid_t uid = geteuid();
    passwd *pw = getpwuid(uid);
    if (pw != nullptr)
    {
      res = pw->pw_name;
    }
#endif

    return std::string(res);
  }

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

  DataExporter::DataExporter(const ExecInfo &info,
                             std::string_view _filename,
                             std::optional<export_metadata_t> user_description)
      : pimpl(new impl(_filename))
  {
    export_metadata_t description = user_description.has_value()
                                        ? *user_description
                                        : "Interesting results";

    metadata["file_version"] = 4;
    metadata["creation_date"] = date_time();
    metadata["author"] = get_user_name();
    metadata["description"] = description;
    metadata["run_id"] = info.run_id;
  }

  void DataExporter::do_link(std::string_view filename,
                             std::string_view groupname)
  {
    pimpl->file.createExternalLink(
        "tes_link", filename.data(), groupname.data());
  }

  void
  DataExporter::write_properties(std::optional<std::string> specific_dataspace,
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

  void DataExporter::prepare_matrix(MultiMatrixDescription description)
  {
    CHECK_PIMPL
    constexpr size_t hdf5_max_compression = 9;

    HighFive::DataSpace dataspace(description.dims, description.max_dims);
    HighFive::DataSetCreateProps props;
    if (description.chunk_dims.has_value())
    {
      props.add(HighFive::Chunking(*description.chunk_dims));
    }
    props.add(HighFive::Shuffle());
    if (description.compression)
    {
      props.add(HighFive::Deflate(hdf5_max_compression));
    }

    if (description.is_integer)
    {
      pimpl->file.createDataSet<uint64_t>(description.name, dataspace, props);
    }
    else
    {
      pimpl->file.createDataSet<double>(description.name, dataspace, props);
    }
  }

  void DataExporter::write_simple(std::string specific_dataspace,
                                  const simple_export_t &value)
  {
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
          pimpl->file.createDataSet<T>(specific_dataspace, val);
        },
        value);
  }

  void DataExporter::write_simple(const export_initial_kv &values)
  {

    for (const auto &kv : values)
    {
      const std::string &attributeName = "initial_parameters/" + kv.first;
      const simple_export_t &value = kv.second;
      write_simple(attributeName, value);
    }
  }

  DataExporter::~DataExporter() = default;
} // namespace DE__WIP
