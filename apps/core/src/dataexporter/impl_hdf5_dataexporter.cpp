#include <stdexcept>
#ifdef USE_HIGHFIVE
#  include <chrono>
#  include <common/common.hpp>
#  include <cstddef>
#  include <ctime>
#  include <dataexporter/data_exporter.hpp>
#  include <highfive/H5DataSpace.hpp>
#  include <highfive/H5File.hpp>
#  include <highfive/H5PropertyList.hpp>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>
#  include <iomanip>
#  include <sstream>
#  include <string>
#  include <string_view>
#  include <type_traits>
#  include <variant>
# include <span> 
#  ifdef __linux__
#    include <pwd.h>
#    include <sys/types.h>
#    include <unistd.h>
#  endif

#  define CHECK_PIMPL                                                                              \
    if (!pimpl)                                                                                    \
    {                                                                                              \
      throw std::runtime_error(__FILE__ ": Unexpected ERROR");                                     \
    }
namespace
{
  std::string date_time()
  {
    // Non vedo l’ora che arrivi c++23-format
    std::stringstream ss;
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
    return ss.str();
  }

  std::string get_user_name()
  {
    std::string_view res = "someone";
#  ifdef __linux__
    uid_t uid = geteuid();
    passwd* pw = getpwuid(uid);
    if (pw != nullptr)
    {
      res = pw->pw_name;
    }
#  endif

    return std::string(res);
  }

  std::size_t get_chunk_size(std::size_t data_length)
  {
    // NOLINTBEGIN
    if (data_length <= 1000)
    {
      return 1;
    }
    else if (data_length <= 1e4)
    {
      return 1024;
    }
    else if (data_length <= 1e6)
    {
      return 8192;
    }
    else
    {
      return 65536;
    }
    // NOLINTEND
  }
} // namespace

namespace Core
{

  class DataExporter::impl
  {
  public:
    explicit impl(std::string_view _filename)
        : filename(_filename), file(HighFive::File(filename, HighFive::File::Truncate))
    {
    }

    impl(const impl&) = delete;
    impl(impl&&) = delete;
    impl& operator=(const impl&) = delete;
    impl& operator=(impl&&) = delete;

    ~impl() = default;

    template <typename T> void write(std::string_view name, T val)
    {
      file.createAttribute<T>(std::string(name), val);
    }

    std::string filename;
    HighFive::File file;
  };

  constexpr size_t hdf5_max_compression = 9;
  DataExporter::DataExporter(const ExecInfo& info,
                             std::string_view _filename,

                             std::optional<export_metadata_t> user_description)
      : pimpl(new impl(_filename))
  {
    export_metadata_t description =
        user_description.has_value() ? *user_description : "Interesting results";

    metadata["file_version"] = 5; //NOLINT
    metadata["creation_date"] = date_time();
    metadata["author"] = get_user_name();
    metadata["description"] = description;
    metadata["run_id"] = info.run_id;
  }

  void DataExporter::do_link(std::string_view filename,
                             std::string_view link_name,
                             std::string_view groupname)
  {
    pimpl->file.createExternalLink(link_name.data(), filename.data(), groupname.data());
  }

  void DataExporter::write_properties([[maybe_unused]] std::optional<std::string> specific_dataspace,
                                      const export_metadata_kv& values)
  {
    // TODO impl specific dataspace as metadata by default.
    // It's not easy to deal with built-in hdf5 metadata from other lib(python/rust),
    // use standard dataset called metadata instead is better
    CHECK_PIMPL

    for (const auto& kv : values)
    {
      const std::string& attributeName = kv.first;
      const export_metadata_t& value = kv.second;
      std::visit(
          [&](const auto& val)
          {
            using T = std::decay_t<decltype(val)>; // Get the actual type T of
                                                   // the variant alternative
            pimpl->write<T>(attributeName, val);
          },
          value);
    }
    pimpl->file.flush();
  }

  void DataExporter::prepare_matrix(MultiMatrixDescription description)
  {
    CHECK_PIMPL

    HighFive::DataSpace dataspace(description.dims, description.max_dims);
    HighFive::DataSetCreateProps props;

    if (description.chunk_dims.has_value())
    {
      auto chunk_dims = description.chunk_dims.value();
      if(chunk_dims.size() != chunk_dims.size())
      {
        throw std::invalid_argument("prepare_matrix(HDF5): container and chunkdimensions and  don´t match");
      }
      props.add(HighFive::Shuffle());
      props.add(HighFive::Chunking(chunk_dims));
    }

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
    descriptors.emplace(description.name, description);
  }

  void DataExporter::write_simple(std::string specific_dataspace, const simple_export_t& value)
  {
    std::visit(
        [&](const auto& val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
          pimpl->file.createDataSet<T>(specific_dataspace, val);
        },
        value);
  }

  void DataExporter::write_simple(const export_initial_kv& values, std::string_view root)
  {

    for (const auto& kv : values)
    {
      std::string path = root.data() + kv.first;
      const simple_export_t& value = kv.second;
      write_simple(path, value);
    }
    pimpl->file.flush();
  }

  void DataExporter::append_array(std::string_view name,
                                  std::span<const double> data,
                                  uint64_t last_size)
  {
    CHECK_PIMPL
    auto& descriptor = descriptors.at(std::string(name));

    auto dataset = pimpl->file.getDataSet(name.data());

    auto dims = descriptor.dims;
    // Following are explicit copy
    auto new_size = dims;
    auto select_start = dims;
    auto select_size = dims;
    new_size[0] = last_size + data.size();
    select_start[0] = (last_size > 0) ? last_size - 1 : 0;
    select_size[0] = data.size();
    dataset.resize(new_size);

    dataset.select(select_start, select_size).write_raw(data.data());
    pimpl->file.flush();
  }

  void
  DataExporter::write_matrix(std::string_view name, std::span<const double> values, bool compress)
  {
    CHECK_PIMPL

    HighFive::DataSetCreateProps ds_props;
    // Minimum chunk size is 1 and leads to error if size<chunk
    // With size=0, nothing will be saved it will not change anything to skip property in this case
    // TODO: if values is 0 early return ?
    if (values.size() > 1)
    {
      // If error occurs, try to debug with fixed chunk of 1
      ds_props.add(HighFive::Chunking(std::min(get_chunk_size(values.size()),values.size())));
      ds_props.add(HighFive::Shuffle());
    }
    const auto data_space = HighFive::DataSpace(values.size());

    if (compress)
    {
      ds_props.add(HighFive::Deflate(hdf5_max_compression));
    }

    auto dataset = pimpl->file.createDataSet<double>(name.data(), data_space, ds_props);
    dataset.write_raw(values.data());
    pimpl->file.flush();
  }

  void DataExporter::write_matrix(std::string_view name,
                                  std::span<const double> values,
                                  size_t n_row,
                                  size_t n_col,
                                  bool compress)
  {
    CHECK_PIMPL
    auto data = Eigen::Map<Eigen::MatrixXd>(
        const_cast<double*>(values.data()), EIGEN_INDEX(n_row), EIGEN_INDEX(n_col));
    HighFive::DataSetCreateProps ds_props;
    ds_props.add(HighFive::Chunking({n_row, n_col}));
    ds_props.add(HighFive::Shuffle());
    if (compress)
    {
      ds_props.add(HighFive::Deflate(hdf5_max_compression));
    }
    pimpl->file.createDataSet(name.data(), data, ds_props);

    pimpl->file.flush();
  }

  void DataExporter::append_matrix(std::string_view name, matrix_variant_t data)
  {
    CHECK_PIMPL
    auto& descriptor = descriptors.at(std::string(name));

    auto dataset = pimpl->file.getDataSet(name.data());

    auto dims = descriptor.dims;
    auto new_size = dims;
    auto select_start = dims;
    auto select_size = dims;
    new_size[0] = (export_counter) + 1;
    select_start[0] = export_counter;
    for (auto i = 1UL; i < select_start.size(); ++i)
    {
      select_start[i] = 0;
    }
    select_size[0] = 1;

    dataset.resize(new_size);

    std::visit(
        [&](auto&& arg)
        {
          using T = std::decay_t<decltype(arg)>;
          if constexpr (std::is_same_v<T, double>)
          {
            dataset.select(select_start, select_size).write_raw(&arg);
          }
          else
          {
            dataset.select(select_start, select_size).write_raw(arg.data());
          }
        },
        data);
    pimpl->file.flush();
  }

  DataExporter::~DataExporter() = default;

} // namespace Core

#endif
