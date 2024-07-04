#include "common/simulation_parameters.hpp"


#include "simulation/simulation.hpp"

#include <data_exporter.hpp>
#include <mc/events.hpp>

#include <Eigen/Core>
#include <chrono>
#include <nl_types.h>
#include <numeric>
#include <string_view>
#include <tuple>

#include <iomanip>

#ifdef USE_HIGHFIVE
#  include "highfive/H5DataSpace.hpp"
#  include "highfive/H5File.hpp"
#  include "highfive/H5PropertyList.hpp"
#  include <Eigen/Dense>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>
#endif

std::string date_time()
{
  std::stringstream ss;
  auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

#ifdef USE_HIGHFIVE

class DataExportHighFive : public DataExporter
{
public:
  DataExportHighFive(const ExecInfo &info,
                     const SimulationParameters &params,
                     std::string_view _filename,
                     std::tuple<size_t, size_t> dim,
                     size_t niter,
                     std::span<size_t> distribution);

  void prepare();
  void append(double t,
              std::span<double> concentration_liquid,
              const std::vector<size_t> &distribution,
              std::span<const double> liquid_volume,
              std::span<const double> volume_gas) override;

  void write_final_particle_data(
      const std::unordered_map<std::string, std::vector<model_properties_t>> &,
      const std::unordered_map<std::string, std::vector<double>> &) override;

protected:
  void write_final_results(ExportData &data,
                           std::span<size_t> distribution) override;

private:
  static void write_attributes(HighFive::File &file, export_metadata_kv &md);

  static void write_initial(HighFive::File &file,
                            const ExecInfo &info,
                            export_initial_kv &md);
};

#endif

DataExporter::DataExporter(const ExecInfo &info,
                           const SimulationParameters &params,
                           std::string_view _filename,
                           std::tuple<size_t, size_t> dim,
                           size_t niter,
                           std::span<size_t> distribution)
    : filename(_filename), n_row(get<0>(dim)), n_col(get<1>(dim)),
      n_iter(niter + 2)

{

  metadata["file_version"] = 3;
  metadata["creation_date"] = date_time();
  metadata["author"] = "someone";
  metadata["description"] = "Interesting results";

  initial_values["number_particles"] = params.user_params.numper_particle;
  initial_values["number_compartment"] = params.n_compartments;
  initial_values["final_time"] = params.user_params.final_time;
  initial_values["particle_distribution"] =
      std::vector<size_t>(distribution.begin(), distribution.end());
  initial_values["delta_time"] = params.d_t;
  initial_values["n_map"] = params.n_different_maps;
  initial_values["t_per_flow_map"] = params.t_per_flow_map;
}

void DataExporter::write_final_results(Simulation::SimulationUnit &simulation,
                                       std::span<size_t> distribution)
{

  size_t n_part = std::accumulate(
      distribution.begin(), distribution.end(), static_cast<size_t>(0));

  ExportData data = {
      n_part,
      simulation.getCliqData(),
      simulation.getCgasData(),
      simulation.mc_unit->ts_events.data(),
      simulation.getDim(),
  };

  write_final_results(data, distribution);
  // Results depending on simulation are exported or copied into properties,
  // clear montecarlo state to save memory
  simulation.clear_mc();
}

std::unique_ptr<DataExporter>
DataExporter::factory(const ExecInfo &info,
                      const SimulationParameters &params,
                      std::string_view _filename,
                      std::tuple<size_t, size_t> dim,
                      size_t niter,
                      std::span<size_t> distribution)
{
#ifdef USE_HIGHFIVE
  return std::make_unique<DataExportHighFive>(
      info, params, _filename, dim, niter, distribution);
#else
  return std::make_unique<DataExporter>(
      info, params, _filename, dim, niter, distribution);
#endif
}



#ifdef USE_HIGHFIVE

void DataExportHighFive::write_attributes(HighFive::File &file,
                                          export_metadata_kv &md)
{

  for (const auto &kv : md)
  {
    const std::string &attributeName = kv.first;
    const export_metadata_t &value = kv.second;
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
                                                 // the variant alternative
          file.createAttribute<T>(attributeName, val);
        },
        value);
  }
}

void DataExportHighFive::write_initial(HighFive::File &file,
                                       const ExecInfo &info,
                                       export_initial_kv &md)
{

  for (const auto &kv : md)
  {
    const std::string &attributeName = "initial_parameters/" + kv.first;
    const export_initial_t &value = kv.second;
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
          // the variant alternative
          file.createDataSet<T>(attributeName, val);
        },
        value);
  }

  file.createDataSet("misc/n_node_thread", info.thread_per_process);
  file.createDataSet("misc/n_rank", info.n_rank);
}

DataExportHighFive::DataExportHighFive(const ExecInfo &info,
                                       const SimulationParameters &params,
                                       std::string_view _filename,
                                       std::tuple<size_t, size_t> dim,
                                       size_t niter,
                                       std::span<size_t> distribution)
    : DataExporter(info, params, _filename, dim, niter, distribution)
{
  HighFive::File file(filename, HighFive::File::Truncate);
  write_initial(file, info, initial_values);
  write_attributes(file, metadata);
  prepare();
}
constexpr size_t hdf5_max_compression = 9;
void DataExportHighFive::prepare()
{
  HighFive::File file(filename, HighFive::File::ReadWrite);

  auto create_dataset = [&file](const std::string &name,
                                const std::vector<size_t> &dims,
                                const std::vector<size_t> &max_dims,
                                const HighFive::DataType &dtype,
                                const std::vector<hsize_t> &chunk_dims,
                                size_t compression = hdf5_max_compression)
  {
    HighFive::DataSpace dataspace(dims, max_dims);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk_dims));
    props.add(HighFive::Shuffle());
    props.add(HighFive::Deflate(compression));
    return file.createDataSet(name, dataspace, dtype, props);
  };

  auto double_type = HighFive::create_datatype<double>();
  auto size_t_type = HighFive::create_datatype<size_t>();
  create_dataset("records/concentration_liquid",
                 {1, n_col, n_row},
                 {n_iter, n_col, n_row},
                 double_type,
                 {1, n_col, n_row});
  create_dataset("records/liquid_volume",
                 {1, n_col},
                 {n_iter, n_col},
                 double_type,
                 {1, n_col});
  create_dataset("records/gas_volume",
                 {1, n_col},
                 {n_iter, n_col},
                 double_type,
                 {1, n_col});
  create_dataset("records/distribution",
                 {1, n_col},
                 {n_iter, n_col},
                 size_t_type,
                 {1, n_col});
  create_dataset("records/time", {1}, {n_iter}, double_type, {1});

  file.flush();
}

void DataExportHighFive::append(double t,
                                std::span<double> concentration_liquid,
                                const std::vector<size_t> &distribution,
                                std::span<const double> liquid_volume,
                                std::span<const double> volume_gas)
{
  try
  {

    HighFive::File file(filename, HighFive::File::ReadWrite);

    auto write_dataset =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            std::span<T> data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data.data());
    };

    auto write_dataset_vector =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            const std::vector<T> &data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data.data());
    };

    auto write_dataset_real =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            const T *data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data);
    };

    write_dataset("records/concentration_liquid",
                  {counter + 1, n_col, n_row},
                  {counter, 0, 0},
                  {1, n_col, n_row},
                  concentration_liquid);
    write_dataset("records/liquid_volume",
                  {counter + 1, n_col},
                  {counter, 0},
                  {1, n_col},
                  liquid_volume);
    write_dataset("records/gas_volume",
                  {counter + 1, n_col},
                  {counter, 0},
                  {1, n_col},
                  volume_gas);
    write_dataset_vector("records/distribution",
                         {counter + 1, n_col},
                         {counter, 0},
                         {1, n_col},
                         distribution);
    write_dataset_real("records/time", {counter + 1}, {counter}, {1}, &t);

    counter++;
  }
  catch (...)
  {
  }
}

void DataExportHighFive::write_final_results(ExportData &data,
                                             std::span<size_t> distribution)

{

  const auto dim = data.dim;
  const auto nr = static_cast<size_t>(get<0>(dim));
  const auto nc = static_cast<size_t>(get<1>(dim));

  HighFive::File file(filename, HighFive::File::ReadWrite);

  file.createDataSet("final_results/number_particles", data.number_particles);

  auto cl =
      Eigen::Map<Eigen::MatrixXd>(data.concentration_liquid.data(), nr, nc);

  file.createDataSet("final_results/concentrations/liquid", cl);

  if (!data.concentration_gas.empty())
  {
    auto cg =
        Eigen::Map<Eigen::MatrixXd>(data.concentration_gas.data(), nr, nc);
    file.createDataSet("final_results/concentrations/gas", cg);
  }

  file.createDataSet(
      "final_results/distribution",
      std::vector<double>(distribution.begin(), distribution.end()));

  file.createDataSet("final_results/events/total_division",
                     data.event->get<MC::EventType::NewParticle>());
  file.createDataSet("final_results/events/total_death",
                     data.event->get<MC::EventType::Death>());

  file.createDataSet<std::string>("biological_model/description ",
                                  std::string("Model description"));
}

void DataExportHighFive::write_final_particle_data(
    const std::unordered_map<std::string, std::vector<model_properties_t>>
        &props,
    const std::unordered_map<std::string, std::vector<double>> &spatial_props)
{
  HighFive::File file(filename, HighFive::File::ReadWrite);
  std::cout << "EXPORTING PARTICLE DATA" << std::endl;
  HighFive::DataSetCreateProps ds_props;
  ds_props.add(HighFive::Chunking(props.size()));
  ds_props.add(HighFive::Shuffle());
  ds_props.add(HighFive::Deflate(hdf5_max_compression));

  for (const auto &[key, values] : props)
  {
    auto v1 = values[0];
    const auto size = values.size();

    std::visit(
        [&file, key = std::move(key), values = std::move(values), size](
            const auto &sample_val)
        {
          using T = std::decay_t<decltype(sample_val)>;
          if constexpr (std::is_same_v<T, double>)
          {
            std::vector<size_t> dim({size});
            std::vector<double> non_zero_values;
            for (const auto &value : values)
            {
              if (std::holds_alternative<double>(value))
              {
                double val = std::get<double>(value);
                if (val != 0.0)
                {
                  non_zero_values.push_back(val);
                }
              }
            }
            auto ds =
                file.createDataSet("biological_model/" + key, non_zero_values);
          }
          else
          {
            return;
          }
        },
        v1);
  }

  HighFive::DataSetCreateProps _ds_props;
  _ds_props.add(HighFive::Chunking(1));
  _ds_props.add(HighFive::Shuffle());
  _ds_props.add(HighFive::Deflate(hdf5_max_compression));
  for (const auto &[key, values] : spatial_props)
  {
    file.createDataSet("biological_model/spatial/" + key, values, _ds_props);
  }
}

#endif
